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

    BufferedTable(std::string name, const table_t &table) :
            table_t(table), m_buffer(name, 1) {
        Table::set_buffer(&m_buffer);
    }

    template<typename ...Args>
    BufferedTable(std::string name, Args... args):
            BufferedTable(name, table_t(args...)) {}

    BufferedTable(const BufferedTable<table_t> &other) :
            BufferedTable(other.m_buffer.m_name, static_cast<const table_t &>(other)) {}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};


#endif //M7_BUFFEREDTABLE_H
