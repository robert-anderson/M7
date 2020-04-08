//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_BITSET_H
#define M7_BITSET_H

#include "Table.h"
#include "BitsetField.h"

class Bitset : public BitsetElement {
    struct BitsetTable : public Table {
        BitsetField field;
        BitsetTable(Bitset* bitset, size_t nbit) :
            Table(),
            field(this, 1, nbit){expand(1);}
    };
    BitsetTable internal_table;
public:
    Bitset(size_t nbit):BitsetElement(nullptr, nullptr), internal_table(this, nbit){
        m_field = &internal_table.field;
        auto element = (*m_field)(0);
        m_begin = element.begin();
    }
};

#endif //M7_BITSET_H
