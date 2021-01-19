//
// Created by rja on 01/11/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include <forward_list>
#include <stack>
#include <src/core/parallel/RankAllocator.h>
#include "Table.h"
#include "src/core/field/TableField.h"

struct LookupResult {
    std::forward_list<size_t>& m_bucket;
    std::forward_list<size_t>::iterator m_prev;
    LookupResult(std::forward_list<size_t>& bucket):
    m_bucket(bucket), m_prev(bucket.before_begin()){}
    operator bool() const {
        return m_prev!=m_bucket.end();
    }
    size_t operator*() const {
        if (!*this) return ~0ul;//throw std::runtime_error("Cannot return stored index, hashmap entry not found!");
        return *std::next(m_prev);
    }
};

template<typename table_t, typename field_t, typename hash_fn=typename field_t::hash_fn>
struct MappedTable : table_t {
    static_assert(std::is_base_of<Table, table_t>::value, "Template arg must be derived from Table");

    using Table::push_back;
    using Table::m_row_size;
    using Table::m_row_dsize;
    using Table::m_bw;
    using Table::clear;
    using Table::dbegin;

    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 2.0;
    static constexpr size_t c_nbucket_min = 100;
    static_assert(std::is_base_of<NdFieldGroup<0ul>, field_t>::value, "Key field must be scalar");

    field_t& m_key_field;
    std::vector<std::forward_list<size_t>> m_buckets;
    std::stack<size_t> m_free_rows;

    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_total_skips = 0.0;
    size_t m_total_lookups = 0.0;
    double m_remap_ratio = c_default_remap_ratio;
    size_t m_remap_nlookup = c_default_remap_nlookup;

    template<typename ...Args>
    MappedTable(field_t &key_field, Args... args) :
            table_t(args...), m_key_field(key_field), m_buckets(c_nbucket_min) {}

    template<typename ...Args>
    MappedTable(size_t nbucket, field_t &key_field, Args... args) :
            table_t(args...), m_key_field(key_field), m_buckets(nbucket) {}

    defs::hash_t hash (typename field_t::view_t key){
        return hash_fn()(key);
    }

    typename field_t::view_t get_key(const size_t& irow){
        return m_key_field(irow);
    }

    LookupResult operator[](const typename field_t::view_t key) {
        m_total_lookups++;
        LookupResult res(m_buckets[hash(key) % m_buckets.size()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next!= res.m_bucket.end(); next++){
            if (get_key(*next)==key) return res;
            res.m_prev++;
            m_total_skips++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    void erase_rows(const defs::inds& irows) override {
        for (auto irow : irows) {
            auto lookup = (*this)[m_key_field(irow)];
            erase(lookup);
        }
    }

    void insert_rows(const Buffer &recv) override {
        const auto nrow = recv.dsize()/m_row_dsize;
        for (size_t irow = 0; irow<nrow; ++irow){
            auto iinsert = get_free_row();
            std::memcpy(dbegin(iinsert), recv.dbegin()+irow*m_row_dsize, m_row_size);
            post_insert(iinsert);
        }
    }

    void erase(LookupResult result){
        clear(*result);
        m_free_rows.push(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }

    size_t get_free_row(){
        if (m_free_rows.empty()) return push_back();
        auto irow = m_free_rows.top();
        m_free_rows.pop();
        return irow;
    }

    size_t insert(const typename field_t::view_t key) {
        ASSERT(!(this->operator[](key)))
        auto irow = get_free_row();
        auto& bucket = m_buckets[hash(key) % m_buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_key_field(irow) = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert(const size_t& irow, std::vector<std::forward_list<size_t>>& buckets) {
        ASSERT(!Table::is_cleared(irow))
        auto key = m_key_field(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[key]);
        auto& bucket = buckets[hash(key) % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t& irow) {
        post_insert(irow, m_buckets);
    }

    void print_map() const {
        size_t nkey = 0ul;
        for (size_t ibucket=0ul; ibucket<m_buckets.size(); ++ibucket){
            auto bucket = m_buckets[ibucket];
            if (bucket.empty()) continue;
            std::cout << "Bucket " << std::to_string(ibucket) << std::endl;
            for (auto irow: bucket){
                std::cout << "\t" << m_key_field(irow).to_string() << " => " << irow << std::endl;
                nkey++;
            }
        }
        std::cout << "Number of keys: " << nkey << std::endl;
    }

    void remap_if_required() {
        if (m_total_lookups<m_remap_nlookup) return;
        if (double(m_total_skips)/double(m_total_lookups)>m_remap_ratio) return;
        // use the same expansion factor as for the Table buffer
        const size_t nbucket_new = m_buckets.size()*m_bw.expansion_factor();
        std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
        for (const auto& old_bucket : m_buckets) {
            for (const auto irow : old_bucket){
                post_insert(irow, new_buckets);
            }
        }
        m_buckets = std::move(new_buckets);
        // reset counters
        m_total_skips = 0ul;
        m_total_lookups = 0ul;
    }

    /**
     * @param nrow_max
     * expected maximum number of rows
     * @param mean_skips_per_bucket
     * the maximum acceptable value for the average number
     * @return
     */
    static size_t nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket) {
        // number of skips is on average half the average number of rows per bin.
        return std::max(c_nbucket_min, size_t(nrow_max/double(2*mean_skips_per_bucket)));
    }

};


#endif //M7_MAPPEDTABLE_H
