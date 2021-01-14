//
// Created by rja on 01/11/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include <forward_list>
#include <stack>
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


template<typename field_t, typename hash_fn=typename field_t::hash_fn>
struct MappedTable : Table {

    static_assert(std::is_base_of<NdFieldGroup<0ul>, field_t>::value, "Key field must be scalar");

    field_t& m_key_field;
    std::vector<std::forward_list<size_t>> m_buckets;
    std::stack<size_t> m_free_rows;

    size_t nbucket() const {
        return m_buckets.size();
    }

    MappedTable(field_t &key_field, size_t nbucket) :
    Table(), m_key_field(key_field), m_buckets(nbucket) {}


    defs::hash_t hash (typename field_t::view_t key){
        return hash_fn()(key);
    }

    typename field_t::view_t get_key(const size_t& irow){
        return m_key_field(irow);
    }

    LookupResult operator[](const typename field_t::view_t key) {
        LookupResult res(m_buckets[hash(key) % nbucket()]);
        auto current = res.m_bucket.before_begin();
        for (auto next = std::next(current); next!= res.m_bucket.end(); next++){
            if (get_key(*next)==key) return res;
            res.m_prev++;
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
        auto& bucket = m_buckets[hash(key) % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_key_field(irow) = key;
        return irow;
    }

    // for situations in which the row has already been copied-to
    void post_insert(const size_t& irow) {
        ASSERT(!Table::is_cleared(irow))
        auto key = m_key_field(irow);
        // row mustn't have already been added
        ASSERT(!(*this)[key]);
        auto& bucket = m_buckets[hash(key) % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
    }

    void print_map() const {
        size_t nkey = 0ul;
        for (size_t ibucket=0ul; ibucket<nbucket(); ++ibucket){
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

};


#endif //M7_MAPPEDTABLE_H
