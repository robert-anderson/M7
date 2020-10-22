//
// Created by rja on 02/10/2020.
//

#ifndef M7_BUFFEREDTABLEARRAY_H
#define M7_BUFFEREDTABLEARRAY_H

#include "Buffer.h"
#include "BufferWindow.h"

template<typename table_t>
class BufferedTableArray {
    Buffer m_buffer;
    std::vector<table_t> m_tables;
    size_t m_nrow_per_table = 0ul;

    void move_tables(Buffer& new_buffer, size_t nrow_per_table, bool expansion=true){
        const auto dsize_per_table = nrow_per_table*defs::ndata_cacheline;
        for (size_t itable=0ul; itable<size(); ++itable){
            // move tables in reverse if expanding the buffer
            auto& table = expansion ? m_tables[size()-1-itable] : m_tables[itable];
            table.move(BufferWindow(new_buffer, dsize_per_table*itable, dsize_per_table));
        }
    }

public:
    template<typename ...Args>
    BufferedTableArray(size_t ntable, Args... args): table_t(args...), m_buffer(){
        table_t::move(BufferWindow(m_buffer));
    }

    size_t size() const {return m_tables.size();}

    void resize(size_t nrow_per_table) {
        Buffer new_buffer(table_t::m_ncacheline, nrow_per_table*size());
        move_tables(new_buffer, nrow_per_table, nrow_per_table>m_nrow_per_table);
        m_nrow_per_table = nrow_per_table;
        m_buffer = std::move(new_buffer);
    }

    void expand(size_t delta_nrow){
        resize(table_t::m_nrow+delta_nrow);
    }

    table_t& operator[](const size_t& itable){
        return m_tables[itable];
    }

    const table_t& operator[](const size_t& itable) const {
        return m_tables[itable];
    }
};



#endif //M7_BUFFEREDTABLEARRAY_H
