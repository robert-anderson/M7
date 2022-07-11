//
// Created by Robert J. Anderson on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "MappedTable.h"

template<typename row_t, bool mapped=false>
class BufferedTable : public std::conditional<mapped, MappedTable<row_t>, Table<row_t>>::type {
    Buffer m_buffer;
public:
    typedef typename std::conditional<mapped, MappedTable<row_t>, Table<row_t>>::type table_t;
    using TableBase::m_bw;

    BufferedTable(str_t name, const table_t& table):
            table_t(table), m_buffer(name, 1) {
        TableBase::set_buffer(&m_buffer);
        ASSERT(static_cast<const Row&>(Table<row_t>::m_row).m_table);
    }

    BufferedTable& operator=(const BufferedTable<row_t, mapped> &other) {
        table_t::operator=(other);
        return *this;
    }

    bool operator==(const BufferedTable<row_t, mapped> &other) const {
        return static_cast<const table_t&>(*this) == static_cast<const table_t&>(other);
    }

    BufferedTable(const BufferedTable<row_t, mapped> &other) : BufferedTable(other.m_buffer.m_name, other){
        *this = other;
        table_t::m_row.restart();
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    double get_expansion_factor() const {
        return m_buffer.m_expansion_factor;
    }
};

#endif //M7_BUFFEREDTABLE_H
