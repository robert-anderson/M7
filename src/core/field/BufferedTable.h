//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "MappedTable.h"

template<typename row_t>
class BufferedTable : public Table<row_t> {
    Buffer m_buffer;
public:
    using TableBase::m_row_dsize;

    BufferedTable(std::string name, const row_t& row):
            Table<row_t>(row), m_buffer(name, 1) {
        TableBase::set_buffer(&m_buffer);
    }
    BufferedTable(const BufferedTable<row_t> &other) :
            BufferedTable(other.m_buffer.m_name, other.m_row){}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};


template<typename row_t>
class BufferedMappedTableZ : public MappedTable<row_t> {
    Buffer m_buffer;
public:
    using TableBase::m_row_dsize;

    BufferedMappedTableZ(std::string name, const row_t& row, size_t nbucket):
            MappedTable<row_t>(row, nbucket), m_buffer(name, 1) {
        TableBase::set_buffer(&m_buffer);
    }
    BufferedMappedTableZ(const MappedTable<row_t> &other) :
            BufferedMappedTableZ(other.m_buffer.m_name, other.nbucket(), other.m_row){}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};

#endif //M7_BUFFEREDTABLE_H
