//
// Created by Robert J. Anderson on 02/09/2021.
//

#include "BosonCoeffs_2.h"

BosonCoeffs_2::BosonCoeffs_2(size_t nmode) : m_nmode(nmode), m_v(utils::integer::npair(nmode*nmode)){}

void BosonCoeffs_2::set(size_t i, size_t j, size_t k, size_t l, defs::ham_t value) {
    m_v.set(index(i, j, k, l), value);
}

defs::ham_t BosonCoeffs_2::get(size_t i, size_t j, size_t k, size_t l) const {
    return m_v[index(i, j, k, l)];
}

defs::ham_t BosonCoeffs_2::phys_element(size_t i, size_t j, size_t k, size_t l) const {
    return get(i, k, j, l);
}
