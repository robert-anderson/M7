//
// Created by rja on 28/08/2021.
//

#include <M7_lib/field/FrmOnvField.h>
#include "FrmBosCoupledCoeffs.h"

size_t FrmBosCoupledCoeffs::index(size_t n, size_t p, size_t q) const {
    if (m_bd.m_nsite==m_nintind_frm){
        // not storing spin-resolved, so convert spin orbital indices to site indices
        p = FrmOnvField::isite(p, m_bd.m_nsite);
        q = FrmOnvField::isite(q, m_bd.m_nsite);
    }
    return n * m_nintind_frm2 + p * m_nintind_frm + q;
}

FrmBosCoupledCoeffs::FrmBosCoupledCoeffs(BasisData bd, bool spin_resolved):
        m_bd(bd), m_nintind_frm(m_bd.m_nsite * (spin_resolved ? 2ul: 1ul)), m_nintind_frm2(m_nintind_frm*m_nintind_frm),
        m_nintind_bos(bd.m_nmode), m_v(m_nintind_bos*m_nintind_frm2){}

void FrmBosCoupledCoeffs::set(const size_t &n, const size_t &p, const size_t &q, const defs::ham_t &value) {
    m_v.set(index(n, p, q), value);
}

const defs::ham_t &FrmBosCoupledCoeffs::get(const size_t &n, const size_t &p, const size_t &q) const {
    return m_v[index(n, p, q)];
}

bool FrmBosCoupledCoeffs::constant_diagonal() const {
    if (m_bd.m_nsite!=m_bd.m_nmode) return false;
    auto v = get(0,0,0);
    for (size_t i=1ul; i<m_bd.m_nspinorb; ++i) {
        auto n = FrmOnvField::isite(i, m_bd.m_nsite);
        if (get(n,i,i)!=v) return false;
    }
    return true;
}
