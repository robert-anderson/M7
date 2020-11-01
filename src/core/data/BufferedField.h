//
// Created by rja on 01/11/2020.
//

#ifndef M7_BUFFEREDFIELD_H
#define M7_BUFFEREDFIELD_H

#include "BufferedTable.h"
#include "NdField.h"

template<typename field_t>
struct BufferedField {
    static_assert(std::is_base_of<FieldBaseX, field_t>::value,
                  "Template arg must be derived from FieldBase");
    struct InternalTable : TableX {
        NdFieldX<field_t, 0ul> m_field;
        template<typename ...Args>
        InternalTable(Args... args):m_field(this, field_t(args...), "Internal field", {}){}
    };
    BufferedTable<InternalTable> m_table;

    template<typename ...Args>
    BufferedField(Args... args) : m_table(args..., NdFormat<0ul>{}){
        m_table.expand(1);
    }

    typename field_t::view_t operator()(){
        return m_table.m_field(0);
    }
};


#endif //M7_BUFFEREDFIELD_H
