//
// Created by RJA on 25/10/2020.
//

#ifndef M7_FLAGFIELD_H
#define M7_FLAGFIELD_H

#include "BitsetField.h"
#include "Table.h"
#include "Flag.h"


#if 0
struct FlagFieldX : public BitsetFieldX<0> {

    using BitsetFieldX<0>::m_table;
    std::vector<Flag*> m_flags;

    size_t get_bit_offset(){
        // called after adding
        ASSERT(m_table->m_fields.back()==this)
        if (m_table->m_fields.size()==1) return 0;
        const auto prev_field = m_table->m_fields[m_table->m_fields.size()-2];
        //if (prev_field->)

//        const auto prev_field = static_cast<const FlagFieldX*>(
//                m_table->m_fields[m_table->m_fields.size()-2]);
//        if(m_table->m_fields[])
    }

    FlagFieldX(TableX* table, std::string description):
    BitsetFieldX<0>(table, {}, 0, description){}

//    template<typename ...Args>
//    typename BitsetFieldX<0>::View::BitView operator()(const size_t &irow) {
//        return View(*this, irow, 0);
//    }
};


#endif //M7_FLAGFIELD_H
#endif //M7_FLAGFIELD_H
