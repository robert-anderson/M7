//
// Created by Robert John Anderson on 2020-03-26.
//

#ifndef SANDBOX2_FLAG_H
#define SANDBOX2_FLAG_H


#include <cstddef>
#include <climits>
#include "src/defs.h"
#include "src/utils.h"
#include "BitsetField.h"


/*
 * FlagFields always have one element, the multidimensionality is implemented
 * at the Flag level
 */

class FlagField;

class FlagElement {
    BitsetElement m_bitset_element;
    const size_t m_ielement;

public:
    FlagElement(BitsetElement bitset_element, const size_t &ielement):
    m_bitset_element(bitset_element), m_ielement(ielement){}

    void operator=(bool v){
        if (v){
            m_bitset_element.set(m_ielement);
        } else {
            m_bitset_element.clr(m_ielement);
        }
    }
    operator bool() { return m_bitset_element.get(m_ielement); }
};

class Flag {
    FlagField *m_field;
    const size_t m_nelement;
    // flag offset is given by (offset of char, bit number in char)
    const defs::pair m_offset;
public:
    Flag(FlagField *field, size_t nelement);

    FlagElement operator()(const size_t &irow, const size_t &isegment, const size_t &ielement);

    friend class FlagField;
};


#endif //SANDBOX2_FLAG_H
