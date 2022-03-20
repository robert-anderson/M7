//
// Created by rja on 02/09/2021.
//

#include "BosonCoeffs_2.h"

BosonCoeffs_2::BosonCoeffs_2(size_t nmode) : m_nmode(nmode), m_v(trig(0, nmode*nmode)){}

void BosonCoeffs_2::set(const size_t &i, const size_t &j, const size_t &k, const size_t &l, const defs::ham_t &value) {
    m_v.set(index(i, j, k, l), value);
}

const defs::ham_t &BosonCoeffs_2::get(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_v[index(i, j, k, l)];
}

const defs::ham_t &BosonCoeffs_2::phys_element(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return get(i, k, j, l);
}
