//
// Created by Robert John Anderson on 2020-01-27.
//

#ifndef M7_DETERMINANT_H
#define M7_DETERMINANT_H

#include <array>
#include <string>
#include "../data/BitfieldNew.h"

class Determinant {
    /*
     * a determinant is defined by a pair of bit strings,
     * one for each spin channel
     */

public:
    std::array<BitfieldNew, 2> m_bitfields;
    Determinant(const size_t &nspatorb);
    Determinant(const size_t &nspatorb, defs::data_t* data1, defs::data_t* data2);
    std::string to_string() const;
    void print() const;
    void zero();
    void set(const size_t &ispat, const size_t &ispin);
    void set(const size_t &i);
    void set(const defs::inds &inds);
    size_t nexcit(const Determinant &other) const;
    bool phase(const Determinant &other) const;
};


#endif //M7_DETERMINANT_H
