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

    Determinant(const size_t &nspatorb, defs::data_t *data1, defs::data_t *data2);

    Determinant(const BitfieldNew &data1, const BitfieldNew &data2);

    Determinant(const Determinant &det) : Determinant(det.nspatorb()) {
        *this = det;
    };

    std::string to_string() const;

    void print() const;

    void zero();

    bool is_zero() const;

    bool is_null() const;

    void set(const size_t &ispat, const size_t &ispin);

    void set(const size_t &i);

    void set(const defs::inds &inds);

    void clr(const size_t &ispat, const size_t &ispin);

    void clr(const size_t &i);

    void clr(const defs::inds &inds);

    size_t nexcit(const Determinant &other) const;

    bool phase(const Determinant &other) const;

    size_t nelec() const;

    size_t nspatorb() const;

    inline int compare(const Determinant &rhs) const;

    Determinant get_excited_det(const defs::inds &removed, const defs::inds &inserted) const;

    Determinant get_excited_det(const size_t &removed, const size_t &inserted) const;

    Determinant get_excited_det(
            const size_t &removed1, const size_t &removed2,
            const size_t &inserted1, const size_t &inserted2
    ) const;

    Determinant &operator=(const Determinant &rhs);

    bool operator==(const Determinant &rhs) const;

    //inline bool operator>(const BitfieldNew& rhs);
    bool partial_phase(const defs::inds &removed, const size_t &nremoved) const;

    bool partial_phase(const Determinant &other) const;
};


#endif //M7_DETERMINANT_H
