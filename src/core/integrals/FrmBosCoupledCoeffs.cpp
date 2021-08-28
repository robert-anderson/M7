//
// Created by rja on 28/08/2021.
//

#include "FrmBosCoupledCoeffs.h"

size_t FrmBosCoupledCoeffs::index(const size_t &n, const size_t &p, const size_t &q) const {
    return n * m_nmode2 + p * m_nmode + q;
}

FrmBosCoupledCoeffs::FrmBosCoupledCoeffs(size_t nmode) : m_nmode(nmode), m_nmode2(nmode*nmode), m_v(m_nmode2*nmode){}

void FrmBosCoupledCoeffs::set(const size_t &n, const size_t &p, const size_t &q, const defs::ham_t &value) {
    m_v.set(index(n, p, q), value);
}

const defs::ham_t &FrmBosCoupledCoeffs::get(const size_t &n, const size_t &p, const size_t &q) const {
    return m_v[index(n, p, q)];
}
