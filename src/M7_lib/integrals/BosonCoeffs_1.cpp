//
// Created by Robert J. Anderson on 02/09/2021.
//

#include "BosonCoeffs_1.h"

BosonCoeffs_1::BosonCoeffs_1(uint_t nmode) : m_nmode(nmode), m_v(m_nmode*m_nmode){}

void BosonCoeffs_1::set(uint_t n, uint_t m, defs::ham_t value) {
    m_v.set(index(n, m), value);
}

defs::ham_t BosonCoeffs_1::get(uint_t n, uint_t m) const {
    return m_v[index(n, m)];
}
