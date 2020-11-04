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
    template<typename ...Args>
    BufferedTable(Args&&... args): table_t(args...), m_buffer(){
        table_t::move(BufferWindow(m_buffer));
    }

    void resize(size_t nrow) {
        Buffer new_buffer(TableX::m_row_dsize, nrow);
        table_t::move(BufferWindow(new_buffer));
        m_buffer = std::move(new_buffer);
    }

    void expand(size_t delta_nrow){
        resize(table_t::m_nrow+delta_nrow);
    }
};


#endif //M7_BUFFEREDTABLE_H
