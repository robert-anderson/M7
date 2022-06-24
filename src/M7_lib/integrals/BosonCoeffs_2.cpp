//
// Created by Robert J. Anderson on 02/09/2021.
//

#include "BosonCoeffs_2.h"

BosonCoeffs_2::BosonCoeffs_2(uint_t nmode) : m_nmode(nmode), m_v(integer::npair(nmode*nmode)){}

void BosonCoeffs_2::set(uint_t i, uint_t j, uint_t k, uint_t l, ham_t value) {
    m_v.set(index(i, j, k, l), value);
}

ham_t BosonCoeffs_2::get(uint_t i, uint_t j, uint_t k, uint_t l) const {
    return m_v[index(i, j, k, l)];
}

ham_t BosonCoeffs_2::phys_element(uint_t i, uint_t j, uint_t k, uint_t l) const {
    return get(i, k, j, l);
}
