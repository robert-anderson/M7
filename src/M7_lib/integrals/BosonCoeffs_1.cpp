//
// Created by rja on 02/09/2021.
//

#include "BosonCoeffs_1.h"

BosonCoeffs_1::BosonCoeffs_1(size_t nmode) : m_nmode(nmode), m_v(m_nmode*m_nmode){}

void BosonCoeffs_1::set(const size_t &n, const size_t &m, const defs::ham_t &value) {
    m_v.set(index(n, m), value);
}

const defs::ham_t &BosonCoeffs_1::get(const size_t &n, const size_t &m) const {
    return m_v[index(n, m)];
}