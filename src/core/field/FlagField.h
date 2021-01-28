//
// Created by rja on 04/11/2020.
//

#ifndef M7_FLAGFIELD_H
#define M7_FLAGFIELD_H

#include "Flag.h"

template<typename set_t>
struct FlagField : NdField<BitsetSpecifier, 0ul>, set_t {
    static_assert(std::is_base_of<FlagSet, set_t>::value, "Template arg must be derived from Flagset");

    template<typename ...Args>
    FlagField(Table* table, std::string description, const set_t& set):
            NdField<BitsetSpecifier, 0ul>(table, {static_cast<const FlagSet&>(set).nbit()}, description),
    set_t(set){
        FlagSet::m_bitset_field = this;
        m_column.m_data.m_details["type"] = "Flagset";
        m_column.m_data.m_details["number of flags"] = std::to_string(static_cast<const FlagSet*>(this)->m_flags.size());
    }
};


#endif //M7_FLAGFIELD_H
