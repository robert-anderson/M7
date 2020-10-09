//
// Created by RJA on 08/10/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <functional>
#include <forward_list>
#include <src/core/hash/Hashing.h>
#include "Table.h"

template<size_t nmapped_field>
class MappedTable : public Table {
    typedef std::function<defs::hash_t(const char* ptr, const size_t& size)> hash_function_t;
    typedef std::pair<FieldBase*, hash_function_t> mapping_pair_t;
    typedef std::array<mapping_pair_t, nmapped_field> mapping_pairs_t;
    const mapping_pairs_t m_mapping_pairs;

    std::vector<std::forward_list<size_t>> m_buckets;

    MappedTable(size_t nbucket, mapping_pairs_t&& mapping_fields):
        Table(), m_mapping_pairs(std::move(mapping_fields)), m_buckets(nbucket){}


    void set_mapping_pairs(typename mapping_pairs_t::iterator iterator) const {}
    template<typename ...Args>
    void set_mapping_pairs(typename mapping_pairs_t::iterator iterator, FieldBase* field, Args...rest) const {
        /*
         * only a field ref is passed, so use default hasher
         */
        *iterator = mapping_pair_t{field, hashing::fnv_hash};
        iterator++;
        set_mapping_pairs(iterator, rest...);
    }

    template<typename ...Args>
    void set_mapping_pairs(typename mapping_pairs_t::iterator iterator, const mapping_pair_t& pair, Args...rest) const {
        /*
         * a hash function is also passed
         */
        iterator++ = pair;
        set_mapping_pairs(iterator, rest...);
    }

    template<typename ...Args>
    mapping_pairs_t mapping_pairs(Args... args){
        static_assert(sizeof...(args)<=nmapped_field, "Too many FieldBase - HashFunction pairs specified");
        static_assert(sizeof...(args)>=nmapped_field, "Too few FieldBase - HashFunction pairs specified");
        mapping_pairs_t res;
        set_mapping_pairs(res.begin(), args...);
        return res;
    }


    inline defs::hash_t single_hash(const size_t& ipair, const char* ptr, const size_t& size) {
        auto &hash_fn = m_mapping_pairs[ipair].second;
        return hash_fn(ptr, size);
    }

    inline defs::hash_t single_hash(const size_t& ipair, const size_t& irow){
        auto &field = m_mapping_pairs[ipair].first;
        return single_hash(ipair, field.begin(irow), field.size());
    }


    defs::hash_t xor_hash(const defs::hash_t& hash) const {
        return hash;
    }

    template<typename ...Args>
    defs::hash_t xor_hash(const defs::hash_t& hash, const char* first, Args... rest) const {
        constexpr size_t ipair = nmapped_field-sizeof...(rest)-1;
        auto &field = m_mapping_pairs[ipair].first;
        return xor_hash(hash^single_hash(ipair, first, field.size()), rest...);
    }

public:
    template<typename ...Args>
    MappedTable(size_t nbucket, Args... args): Table(),
    m_mapping_pairs(mapping_pairs(args...)),
    m_buckets(nbucket){}

    /*
    defs::hash_t hash(const size_t irow){
        defs::hash_t res = 0ul;
        for (size_t i=0ul; i<nmapped_field; ++i) res|=single_hash(i, irow);
        return res;
    }
     */

/*

    template<typename ...Args>
    defs::hash_t hash(Args... ptrs){
        static_assert(sizeof...(ptrs)<=nmapped_field, "Too many field view char pointers specified");
        static_assert(sizeof...(ptrs)>=nmapped_field, "Too few field view char pointers specified");
        return xor_hash(0ul, ptrs...);
    }


    template<typename ...Args>
    defs::hash_t hash(Args... ptrs){
        static_assert(sizeof...(ptrs)<=nmapped_field, "Too many field view char pointers specified");
        static_assert(sizeof...(ptrs)>=nmapped_field, "Too few field view char pointers specified");
        return xor_hash(0ul, ptrs...);
    }
*/
};


#endif //M7_MAPPEDTABLE_H
