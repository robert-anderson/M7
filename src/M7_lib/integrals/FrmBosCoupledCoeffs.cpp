//
// Created by rja on 28/08/2021.
//

#include <M7_lib/field/FrmOnvField.h>
#include "FrmBosCoupledCoeffs.h"

size_t FrmBosCoupledCoeffs::index(size_t n, size_t p, size_t q) const {
    if (m_bd.m_frm.m_nsite==m_nintind_frm){
        // not storing spin-resolved, so convert spin orbital indices to site indices
        p = FrmOnvField::isite(p, m_bd.m_frm.m_nsite);
        q = FrmOnvField::isite(q, m_bd.m_frm.m_nsite);
    }
    return n * m_nintind_frm2 + p * m_nintind_frm + q;
}

FrmBosCoupledCoeffs::FrmBosCoupledCoeffs(const BasisData& bd, bool spin_resolved):
        m_bd(bd), m_nintind_frm(spin_resolved ? bd.m_frm.m_nspinorb : bd.m_frm.m_nsite),
        m_nintind_frm2(m_nintind_frm*m_nintind_frm), m_nintind_bos(bd.m_bos.m_nmode),
        m_v(m_nintind_bos*m_nintind_frm2){}

void FrmBosCoupledCoeffs::set(const size_t &n, const size_t &p, const size_t &q, const defs::ham_t &value) {
    m_v.set(index(n, p, q), value);
}

const defs::ham_t &FrmBosCoupledCoeffs::get(const size_t &n, const size_t &p, const size_t &q) const {
    return m_v[index(n, p, q)];
}