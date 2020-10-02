//
// Created by rja on 02/10/2020.
//

#ifndef M7_TABLE_H
#define M7_TABLE_H



#include <vector>
#include <climits>
#include <cstring>
#include "BufferWindow.h"

class FieldBase;

struct Table {
    BufferWindow m_bw;
    // length of each row in units of the cacheline length
    size_t m_ncacheline = 0ul;
    size_t m_row_size = 0ul;
    size_t m_tight_row_size = 0ul;
    std::vector<FieldBase*> m_fields;
    size_t m_nrow = 0ul;
    /*
     * "high water mark" is result of the next call to push_back
     */
    size_t m_hwm = 0ul;


    Table(BufferWindow buffer);

    Table();

    void move(BufferWindow new_bw);

    size_t push_back();

    void add_field(FieldBase* field);

    char* row_ptr(const size_t& irow) const;

    void clear();

    void clear_row(const size_t& irow);

    virtual std::string to_string(std::string delimiter="") const;
    virtual std::string to_bit_string(std::string delimiter="") const;
    void print() const;
};



#endif //M7_TABLE_H
