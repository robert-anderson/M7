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
#include "src/core/field/Field.h"
#include "src/core/table/BufferedTable.h"

struct LookupResult {
    std::forward_list<size_t> &m_bucket;
    std::forward_list<size_t>::iterator m_prev;

    LookupResult(std::forward_list<size_t> &bucket) :
            m_bucket(bucket), m_prev(bucket.before_begin()) {}

    operator bool() const {
        return m_prev != m_bucket.end();
    }

    size_t operator*() const {
        if (!*this) return ~0ul;//throw std::runtime_error("Cannot return stored index, hashmap entry not found!");
        return *std::next(m_prev);
    }
};

struct MappedTableBase {
    static constexpr double c_default_remap_ratio = 2.0;
    static constexpr size_t c_default_remap_nlookup = 2.0;
    static constexpr size_t c_nbucket_min = 100;

    std::vector<std::forward_list<size_t>> m_buckets;

    /*
     * skips are counted to avoid having to loop over all buckets to count entries
     */
    size_t m_total_skips = 0.0;
    size_t m_total_lookups = 0.0;
    double m_remap_ratio = c_default_remap_ratio;
    size_t m_remap_nlookup = c_default_remap_nlookup;

protected:
    MappedTableBase(size_t nbucket) : m_buckets(nbucket) {}

public:
    /**
     * @param nrow_max
     * expected maximum number of rows
     * @param mean_skips_per_bucket
     * the maximum acceptable value for the average number
     * @return
     */
    static size_t nbucket_guess(size_t nrow_max, size_t mean_skips_per_bucket) {
        // number of skips is on average half the average number of rows per bin.
        return std::max(c_nbucket_min, size_t(nrow_max / double(2 * mean_skips_per_bucket)));
    }
};

template<typename table_tt, typename key_field_tt, typename hash_fnt=typename key_field_tt::hash_fn>
struct MappedTable : MappedTableBase, table_tt {
    static_assert(std::is_base_of<Table, table_tt>::value, "Template arg must be derived from Table");
    typedef table_tt table_t;
    typedef key_field_tt key_field_t;
    typedef hash_fnt hash_fn;
    using Table::push_back;
    using Table::get_free_row;
    using Table::m_row_size;
    using Table::m_row_dsize;
    using Table::clear;
    using Table::dbegin;
    using Table::m_nrow;
    using Table::m_hwm;

    static_assert(std::is_base_of<NdFieldGroup<0ul>, key_field_t>::value, "Key field must be scalar");

    key_field_t &m_key_field;

private:
    /*
     * key field validation:
     * we assert that m_key_field is a ref to a member of table, i.e. the
     * memory offset between the table and the field is less than the size
     * of table_t
     */
    size_t key_field_offset() const {
        auto kf_offset = std::distance((const char *) this, (const char *) &m_key_field);
        MPI_REQUIRE_ALL(kf_offset > 0 && (size_t) kf_offset <= sizeof(table_t),
                        "Field chosen as the key field must be a member of the table");
        return kf_offset;
    }

    key_field_t &get_key_field_from_offset(size_t nbyte) const {
        return *(key_field_t *) ((char *) (this) + nbyte);
    }

public:
    template<typename ...Args>
    MappedTable(size_t nbucket, key_field_t &key_field, const table_t &table) :
            MappedTableBase(nbucket), table_t(table), m_key_field(key_field) {
        key_field_offset();
    }

    template<typename ...Args>
    MappedTable(key_field_t &key_field, const table_t &table) :
            MappedTable(c_nbucket_min, key_field, table) {}

    MappedTable(const MappedTable &other) :
            MappedTable(get_key_field_from_offset(other.key_field_offset()), other) {}

    defs::hash_t hash(typename key_field_t::view_t key) {
        return hash_fn()(key);
    }

    typename key_field_t::view_t get_key(const size_t &irow) {
        return m_key_field(irow);
    }

    LookupResult operator[](const typename key_field_t::view_t& key) {
        m_total_lookups++;
        LookupResult res(m_buckets[hash(key) % m_buckets.size()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next != res.m_bucket.end(); next++) {
            if (get_key(*next) == key) return res;
            res.m_prev++;
            m_total_skips++;
        }
        res.m_prev = res.m_bucket.end();
        return res;
    }

    void erase_rows(const defs::inds &irows) override {
        for (auto irow : irows) {
            auto lookup = (*this)[m_key_field(irow)];
            erase(lookup);
        }
    }

    void erase(LookupResult result) {
        clear(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }


    size_t insert(const typename key_field_t::view_t key) {
        ASSERT(!(this->operator[](key)))
        auto irow = get_free_row();
        auto &bucket = m_buckets[hash(key) % m_buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_key_field(irow) = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert_buckets(const size_t &irow, std::vector<std::forward_list<size_t>> &buckets) {
        ASSERT(!Table::is_cleared(irow))
        auto key = m_key_field(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[key]);
        auto &bucket = buckets[hash(key) % buckets.size()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void post_insert(const size_t &irow) override {
        post_insert_buckets(irow, m_buckets);
    }

    void print_map() const {
        size_t nkey = 0ul;
        for (size_t ibucket = 0ul; ibucket < m_buckets.size(); ++ibucket) {
            auto bucket = m_buckets[ibucket];
            if (bucket.empty()) continue;
            std::cout << "Bucket " << std::to_string(ibucket) << std::endl;
            for (auto irow: bucket) {
                std::cout << "\t" << m_key_field(irow).to_string() << " => " << irow << std::endl;
                nkey++;
            }
        }
        std::cout << "Number of keys: " << nkey << std::endl;
    }

    void remap_if_required() {
        if (m_total_lookups < m_remap_nlookup) return;
        if (double(m_total_skips) / double(m_total_lookups) > m_remap_ratio) return;
        // use the same expansion factor as for the Table buffer
        const size_t nbucket_new = m_buckets.size() * Table::m_bw.expansion_factor();
        std::vector<std::forward_list<size_t>> new_buckets(nbucket_new);
        for (const auto &old_bucket : m_buckets) {
            for (const auto irow : old_bucket) {
                post_insert(irow, new_buckets);
            }
        }
        m_buckets = std::move(new_buckets);
        // reset counters
        m_total_skips = 0ul;
        m_total_lookups = 0ul;
    }
};


#endif //M7_MAPPEDTABLE_H
