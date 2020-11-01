//
// Created by rja on 01/11/2020.
//

#ifndef M7_FLAGFIELD_H
#define M7_FLAGFIELD_H

#include "BitsetField.h"

struct FlagField;

struct Flag {
    size_t m_nbit;
    size_t m_offset;
    Flag(FlagField* field, size_t nbit);
};

template<size_t nind>
struct NdFlag {

};

struct FlagField {

};


#endif //M7_FLAGFIELD_H
