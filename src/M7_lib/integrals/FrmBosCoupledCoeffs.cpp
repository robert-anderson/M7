//
// Created by Robert J. Anderson on 28/08/2021.
//

#include <M7_lib/field/FrmOnvField.h>
#include "FrmBosCoupledCoeffs.h"

uint_t FrmBosCoupledCoeffs::index(uint_t n, uint_t p, uint_t q) const {
    if (m_sizes.m_frm.m_nsite == m_ncoeff_ind_frm){
        // not storing spin-resolved, so convert spin orbital indices to site indices
        p = m_sizes.m_frm.isite(p);
        q = m_sizes.m_frm.isite(q);
    }
    return n * m_ncoeff_ind_frm2 + p * m_ncoeff_ind_frm + q;
}

FrmBosCoupledCoeffs::FrmBosCoupledCoeffs(sys::Size sizes, bool spin_resolved):
        m_sizes(sizes), m_ncoeff_ind_frm(m_sizes.m_frm.ncoeff_ind(spin_resolved)),
        m_ncoeff_ind_frm2(m_ncoeff_ind_frm * m_ncoeff_ind_frm), m_ncoeff_ind_bos(m_sizes.m_bos.m_nmode),
        m_v(m_ncoeff_ind_bos * m_ncoeff_ind_frm2){}

void FrmBosCoupledCoeffs::set(uint_t n, uint_t p, uint_t q, ham_t value) {
    m_v.set(index(n, p, q), value);
}

ham_t FrmBosCoupledCoeffs::get(uint_t n, uint_t p, uint_t q) const {
    return m_v[index(n, p, q)];
}