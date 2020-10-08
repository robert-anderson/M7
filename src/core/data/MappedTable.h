//
// Created by RJA on 08/10/2020.
//

#ifndef M7_MAPPEDTABLE_H
#define M7_MAPPEDTABLE_H

#include <functional>
#include <forward_list>
#include "Table.h"

template<size_t nmapped_field>
class MappedTable : public Table {
    typedef std::function<defs::hash_t(const char* ptr, const size_t& size)> hash_function_t;
    typedef std::pair<FieldBase&, hash_function_t> mapping_pair_t;
    typedef std::array<mapping_pair_t, nmapped_field> mapping_pairs_t;
    const mapping_pairs_t m_mapping_pairs;

    std::vector<std::forward_list<size_t>> m_buckets;

    MappedTable(size_t nbucket, mapping_pairs_t&& mapped_fields):
        Table(), m_mapping_pairs(std::move(mapped_fields)), m_buckets(nbucket){}

    template<typename ...Args>
    void variadic_to_array_element(
            typename mapping_pairs_t::iterator iterator, const mapping_pair_t& mapping_pair, Args... rest){
        iterator++=mapping_pair;
        variadic_to_array_element(iterator, rest...);
    }

    template<typename ...Args>
    mapping_pairs_t variadic_to_array(Args... args){
        mapping_pairs_t mapping_pairs;
        variadic_to_array_element(mapping_pairs.begin(), args...);
    }

    inline defs::hash_t single_hash(const size_t& ipair, const char* ptr, const size_t& size) {
        auto &hash_fn = m_mapping_pairs[ipair].second;
        return hash_fn(ptr, size);
    }

    inline defs::hash_t single_hash(const size_t& ipair, const size_t& irow){
        auto &field = m_mapping_pairs[ipair].first;
        return single_hash(ipair, field.begin(irow), field.size());
    }

public:
    template<typename ...Args>
    MappedTable(size_t nbucket, Args... args):MappedTable(nbucket, variadic_to_array(args...)){
        static_assert(sizeof...(args)<=nmapped_field, "Too many FieldBase - HashFunction pairs specified");
        static_assert(sizeof...(args)>=nmapped_field, "Too few FieldBase - HashFunction pairs specified");
    }

    /*
    defs::hash_t hash(const size_t irow){
        defs::hash_t res = 0ul;
        for (size_t i=0ul; i<nmapped_field; ++i) res|=single_hash(i, irow);
        return res;
    }
     */

    template<typename ...Args>
    void xor_hash(defs::hash_t& hash, const char* first, Args... rest) const {
        hash^=single_hash(sizeof...(rest))
        xor_hash()
    }


    template<typename ...Args>
    defs::hash_t hash(Args... rest){
        defs::hash_t res = 0ul;
        res = hash_one
        for (size_t i=0ul; i<nmapped_field; ++i) res|=single_hash(i, irow);
        return res;
    }


};


#endif //M7_MAPPEDTABLE_H
