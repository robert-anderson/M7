//
// Created by rja on 02/09/2021.
//

#include "BosonCoeffs.h"

BosonCoeffs::BosonCoeffs(size_t nmode) : m_nmode(nmode), m_v(m_nmode*m_nmode){}

void BosonCoeffs::set(const size_t &n, const size_t &m, const defs::ham_t &value) {
    m_v.set(index(n, m), value);
}

const defs::ham_t &BosonCoeffs::get(const size_t &n, const size_t &m) const {
    return m_v[index(n, m)];
}
