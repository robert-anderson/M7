//
// Created by rja on 01/11/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <unordered_map>
#include <forward_list>
#include <stack>
#include "Table.h"
#include "NdField.h"

struct LookupResult {
    std::forward_list<size_t>& m_bucket;
    std::forward_list<size_t>::iterator m_prev;
    LookupResult(std::forward_list<size_t>& bucket):
    m_bucket(bucket), m_prev(bucket.before_begin()){}
    operator bool() const {
        return m_prev!=m_bucket.end();
    }
    const size_t& operator*() const {
        if (!*this) throw std::runtime_error("Cannot return stored index, hashmap entry not found!");
        return *std::next(m_prev);
    }
};


template<typename field_t, typename hash_fn=typename field_t::hash_fn>
struct MappedTable : TableX {

    static_assert(std::is_base_of<NdFieldGroup<0ul>, field_t>::value, "Key field must be scalar");

    field_t& m_key_field;
    std::vector<std::forward_list<size_t>> m_buckets;
    std::stack<size_t> m_free_rows;

    size_t nbucket() const {
        return m_buckets.size();
    }

    MappedTable(field_t &key_field, size_t nbucket) :
    TableX(), m_key_field(key_field), m_buckets(nbucket) {}


    defs::hash_t hash (typename field_t::view_t key){
        return hash_fn()(m_key_field, key);
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

    void erase(LookupResult result){
        clear_row(*result);
        m_free_rows.push(*result);
        result.m_bucket.erase_after(result.m_prev);
        // put into "not found" state:
        result.m_prev = result.m_bucket.end();
    }

    size_t insert(const typename field_t::view_t key) {
        ASSERT(!(this->operator[](key)))
        size_t irow;
        if (m_free_rows.empty()) irow=push_back();
        else {
            irow = m_free_rows.top(); m_free_rows.pop();
        }
        auto& bucket = m_buckets[hash(key) % nbucket()];
        bucket.insert_after(bucket.before_begin(), irow);
        m_key_field(irow) = key;
        return irow;
    }

};


#endif //M7_MAPPEDTABLE_H
