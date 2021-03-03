//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "MappedTable.h"

template<typename row_t, bool mapped=false>
class BufferedTable : public std::conditional<mapped, MappedTable<row_t>, Table<row_t>>::type {
    Buffer m_buffer;
public:
    typedef typename std::conditional<mapped, MappedTable<row_t>, Table<row_t>>::type table_t;
    using TableBase::m_row_dsize;

    BufferedTable(std::string name, const table_t& table):
            table_t(table), m_buffer(name, 1) {
        TableBase::set_buffer(&m_buffer);
        ASSERT(static_cast<const Row&>(Table<row_t>::m_row).m_table_bw);
    }

    BufferedTable(const BufferedTable<row_t> &other) :
            BufferedTable(other.m_buffer.m_name, other){}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};

#endif //M7_BUFFEREDTABLE_H
