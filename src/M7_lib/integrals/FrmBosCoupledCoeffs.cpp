//
// Created by rja on 28/08/2021.
//

#include <M7_lib/field/FrmOnvField.h>
#include "FrmBosCoupledCoeffs.h"

size_t FrmBosCoupledCoeffs::index(size_t n, size_t p, size_t q) const {
    if (m_extents.m_sites == m_ncoeff_ind_frm){
        // not storing spin-resolved, so convert spin orbital indices to site indices
        p = m_extents.m_sites.isite(p);
        q = m_extents.m_sites.isite(q);
    }
    return n * m_ncoeff_ind_frm2 + p * m_ncoeff_ind_frm + q;
}

FrmBosCoupledCoeffs::FrmBosCoupledCoeffs(BasisExtents extents, bool spin_resolved):
        m_extents(extents), m_ncoeff_ind_frm(m_extents.m_sites.ncoeff_ind(spin_resolved)),
        m_ncoeff_ind_frm2(m_ncoeff_ind_frm * m_ncoeff_ind_frm), m_ncoeff_ind_bos(m_extents.m_nmode),
        m_v(m_ncoeff_ind_bos * m_ncoeff_ind_frm2){}

void FrmBosCoupledCoeffs::set(const size_t &n, const size_t &p, const size_t &q, const defs::ham_t &value) {
    m_v.set(index(n, p, q), value);
}

const defs::ham_t &FrmBosCoupledCoeffs::get(const size_t &n, const size_t &p, const size_t &q) const {
    return m_v[index(n, p, q)];
}