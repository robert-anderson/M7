//
// Created by Robert J. Anderson on 12/9/21.
//

#include "HolsteinLadderHam.h"

ham_t HolsteinLadderHam::get_coeff_1110(uint_t imode, uint_t i, uint_t j) const {
    if (imode != m_basis.m_frm.isite(i)) return 0;
    if (imode != m_basis.m_frm.isite(j)) return 0;
    return m_g;
}

ham_t HolsteinLadderHam::get_coeff_1101(uint_t imode, uint_t i, uint_t j) const {
    return get_coeff_1110(imode, j, i);
}

ham_t HolsteinLadderHam::get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    const auto imode = conn.m_bos.m_cre[0].m_imode;
    const auto nocc_frm = onv.m_frm.site_nocc(imode);
    const auto occ_fac = onv.m_bos.occ_fac(conn.m_bos);
    return m_g*nocc_frm*occ_fac;
}

ham_t HolsteinLadderHam::get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    const auto imode = conn.m_bos.m_ann[0].m_imode;
    const auto nocc_frm = onv.m_frm.site_nocc(imode);
    const auto occ_fac = onv.m_bos.occ_fac(conn.m_bos);
    return m_g*nocc_frm*occ_fac;
}