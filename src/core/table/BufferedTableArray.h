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
        auto new_buffer_dsize = new_buffer.dsize();
        auto new_bw_dsize = new_buffer_dsize/ntable();
        for (size_t itable=0ul; itable<ntable(); ++itable){
            // move tables in reverse if expanding the buffer
            size_t ind = expansion ? ntable()-1-itable : itable;
            auto& table = m_tables[ind];
            table.move(BufferWindow(new_buffer, new_bw_dsize*ind, new_bw_dsize));
        }
    }

public:

    const size_t& nrow_per_table() const {
        return m_nrow_per_table;
    }

    defs::data_t* ptr() {
        return m_buffer.ptr();
    }

    const defs::data_t* ptr() const {
        return m_buffer.ptr();
    }

    size_t buffer_dsize() const {
        return m_buffer.dsize();
    }

    size_t bw_dsize() const {
        return (*this)[0].bw_dsize();
    }

    template<typename ...Args>
    BufferedTableArray(size_t ntable, Args... args): m_buffer(){
        m_tables.reserve(ntable);
        for (size_t itable=0ul; itable<ntable; ++itable) {
            m_tables.emplace_back(args...);
            m_tables.back().move(BufferWindow(m_buffer, 0, 0));
        }
    }

    size_t ntable() const {return m_tables.size();}

    void resize(size_t nrow_per_table) {
        Buffer new_buffer(row_dsize(), nrow_per_table*ntable());
        move_tables(new_buffer, nrow_per_table, nrow_per_table>m_nrow_per_table);
        ASSERT(m_tables[0].m_nrow==nrow_per_table)
        m_buffer = std::move(new_buffer);
    }

    void expand(size_t delta_nrow){
        resize(m_nrow_per_table+delta_nrow);
    }

    table_t& operator[](const size_t& itable){
        return m_tables[itable];
    }

    const table_t& operator[](const size_t& itable) const {
        return m_tables[itable];
    }

    defs::inds hwms() const {
        defs::inds res(ntable());
        for (size_t i=0ul; i<ntable(); ++i) res[i] = (*this)[i].m_hwm;
        return res;
    }

    defs::inds displs() const {
        defs::inds res(ntable());
        for (size_t i=0ul; i<ntable(); ++i) res[i] = buffer_dsize()*i;
        return res;
    }

    void clear() {
        for (auto& table:m_tables) {
            static_cast<TableX&>(table).clear();
            static_cast<TableX&>(table).m_hwm = 0;
        }
    }
};


#endif //M7_BUFFEREDTABLEARRAY_H
