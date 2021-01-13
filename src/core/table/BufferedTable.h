//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFEREDTABLE_H
#define M7_BUFFEREDTABLE_H

#include "BufferWindow.h"
#include "Table.h"

template<typename table_t>
class BufferedTable : public table_t {
    Buffer m_buffer;
public:
    using Table::m_row_dsize;

    template<typename ...Args>
    BufferedTable(std::string name, Args&&... args): table_t(args...),
    m_buffer(name, 0, 0){
        table_t::move(BufferWindow(m_buffer));
    }

    size_t buffer_dsize() const {
        return m_buffer.dsize();
    }

    void resize(size_t nrow) {
        Buffer new_buffer(m_buffer.m_name, Table::m_row_dsize, nrow);
        table_t::move(BufferWindow(new_buffer));
        m_buffer = std::move(new_buffer);
        ASSERT(buffer_dsize() == Table::m_row_dsize*nrow)
    }

    void expand(size_t delta_nrow){
        resize(table_t::m_nrow+delta_nrow);
    }

    void expand_by_factor(double factor){
        expand(std::ceil(factor*table_t::m_nrow));
    }
};


#endif //M7_BUFFEREDTABLE_H
