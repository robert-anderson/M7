//
// Created by RJA on 26/10/2020.
//

#ifndef M7_BUFFEREDTABLEARRAY_H
#define M7_BUFFEREDTABLEARRAY_H



#include "Table.h"

template<typename table_t>
class BufferedTableArray {
    Buffer m_buffer;
    std::vector<table_t> m_tables;
    size_t m_nrow_per_table = 0ul;

    const size_t& row_dsize() const {
        return static_cast<const TableX&>(m_tables[0]).m_row_dsize;
    }


    void move_tables(Buffer& new_buffer, size_t nrow_per_table, bool expansion=true){
        m_nrow_per_table = nrow_per_table;
        for (size_t itable=0ul; itable<size(); ++itable){
            // move tables in reverse if expanding the buffer
            auto& table = expansion ? m_tables[size()-1-itable] : m_tables[itable];
            table.move(BufferWindow(new_buffer, buffer_dsize()*itable, buffer_dsize()));
        }
    }

public:

    defs::data_t* ptr() {
        return m_buffer.ptr();
    }

    const defs::data_t* ptr() const {
        return m_buffer.ptr();
    }

    size_t buffer_dsize() {
        return m_nrow_per_table*row_dsize();
    }

    template<typename ...Args>
    BufferedTableArray(size_t ntable, Args... args): table_t(args...), m_buffer(){
        table_t::move(BufferWindow(m_buffer));
    }

    size_t size() const {return m_tables.size();}

    void resize(size_t nrow_per_table) {
        Buffer new_buffer(row_dsize(), nrow_per_table);
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

    defs::inds hwms() const {
        defs::inds res(size());
        for (size_t i=0ul; i<size(); ++i) res[i] = (*this)[i].m_hwm;
        return res;
    }

    defs::inds displs() const {
        defs::inds res(size());
        for (size_t i=0ul; i<size(); ++i) res[i] = buffer_dsize()*i;
        return res;
    }

    void clear() {
        for (auto& table:m_tables) {
            auto t = static_cast<TableX&>(table);
            t.clear();
            t.m_hwm = 0;
        }
    }
};


#endif //M7_BUFFEREDTABLEARRAY_H
