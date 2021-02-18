//
// Created by rja on 10/02/2021.
//

#ifndef M7_BUFFEREDTABLEZ_H
#define M7_BUFFEREDTABLEZ_H

#include "MappedTableZ.h"

template<typename row_t>
class BufferedTableZ : public TableZ<row_t> {
    Buffer m_buffer;
public:
    using TableBaseZ::m_row_dsize;

    BufferedTableZ(std::string name, const row_t& row):
    TableZ<row_t>(row), m_buffer(name, 1) {
        TableBaseZ::set_buffer(&m_buffer);
    }
    BufferedTableZ(const BufferedTableZ<row_t> &other) :
            BufferedTableZ(other.m_buffer.m_name, other.m_row){}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};


template<typename row_t>
class BufferedMappedTableZ : public MappedTableZ<row_t> {
    Buffer m_buffer;
public:
    using TableBaseZ::m_row_dsize;

    BufferedMappedTableZ(std::string name, const row_t& row, size_t nbucket):
            MappedTableZ<row_t>(row, nbucket), m_buffer(name, 1) {
        TableBaseZ::set_buffer(&m_buffer);
    }
    BufferedMappedTableZ(const MappedTableZ<row_t> &other) :
            BufferedMappedTableZ(other.m_buffer.m_name, other.nbucket(), other.m_row){}

    void set_expansion_factor(double f) {
        m_buffer.m_expansion_factor = f;
    }
};

#endif //M7_BUFFEREDTABLEZ_H
