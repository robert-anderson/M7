//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "Table.h"

template<typename table_t>
class BufferedTable : public table_t {
    static_assert(std::is_base_of<Table, table_t>::value, "Template arg must be derived from Table");
    Buffer m_buffer;
public:
    using Table::m_row_dsize;

    template<typename ...Args>
    BufferedTable(std::string name, Args&&... args): table_t(args...),
    m_buffer(name, 1, 0){
        Table::set_buffer(&m_buffer);
    }

    void set_expansion_factor(double f){
        m_buffer.m_expansion_factor = f;
    }

};


#endif //M7_BUFFEREDTABLE_H
