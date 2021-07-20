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
    using TableBase::m_row_size;

    BufferedTable(std::string name, const table_t& table):
            table_t(table), m_buffer(name, 1) {
        TableBase::set_buffer(&m_buffer);
        ASSERT(static_cast<const Row&>(Table<row_t>::m_row).m_table);
    }

    BufferedTable& operator=(const BufferedTable<row_t, mapped> &other) {
        if (other.m_buffer.size()) {
            m_buffer.resize(other.m_buffer.size());
            Table<row_t>::clear();
            Table<row_t>::push_back(other.m_nrow);
            Table<row_t>::m_bw = other.m_bw;
        }
        return *this;
    }

    BufferedTable(const BufferedTable<row_t, mapped> &other) :
            BufferedTable(other.m_buffer.m_name, other){
        *this = other;
    }

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }

    double get_expansion_factor() const {
        return m_buffer.m_expansion_factor;
    }
};

#endif //M7_BUFFEREDTABLE_H
