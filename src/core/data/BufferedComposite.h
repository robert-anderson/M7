//
// Created by rja on 01/11/2020.
//

#ifndef M7_BUFFEREDCOMPOSITE_H
#define M7_BUFFEREDCOMPOSITE_H

#include "CompositeField.h"
#include "BufferedTable.h"
#include "src/core/nd/NdFormat.h"

template<typename composite_t>
struct BufferedComposite {
    static_assert(std::is_base_of<CompositeField, composite_t>::value,
                  "Template arg must be derived from CompositeField");
    struct InternalTable : TableX {
        composite_t m_composite;
        template<typename ...Args>
        InternalTable(Args... args): m_composite(this, args...){}
    };
    BufferedTable<InternalTable> m_table;

    template<typename ...Args>
    BufferedComposite(Args... args): m_table(args...) {
        m_table.expand(1);
        m_table.push_back();
    }

    typename composite_t::view_t operator()(){
        return m_table.m_composite(0);
    }
};


#endif //M7_BUFFEREDCOMPOSITE_H
