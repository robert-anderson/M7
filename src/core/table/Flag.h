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
class Flag;

class FlagElement {
    BitsetElement m_bitset_element;
    Flag* m_flag;

public:
    FlagElement(Flag *flag, BitsetElement bitset_element);

    void operator=(bool v);

    operator bool();
};

class Flag {
    FlagField *m_field;
    const size_t m_nelement;
    // flag offset is given by (offset of char, bit number in char)
    const defs::pair m_offset;
    const std::string m_description;
public:
    Flag(FlagField *field, size_t nelement, const std::string& description="");

    FlagElement operator()(const size_t &irow, const size_t &isegment=0);

    const defs::pair &offset() const {return m_offset;}

    friend class FlagField;
};


#endif //SANDBOX2_FLAG_H
