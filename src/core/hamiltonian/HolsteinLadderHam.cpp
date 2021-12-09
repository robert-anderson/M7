//
// Created by anderson on 12/9/21.
//

#include "HolsteinLadderHam.h"

defs::ham_t HolsteinLadderHam::get_coeff_0010(const size_t &imode) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_coeff_0001(const size_t &imode) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_coeff_1110(const size_t &imode, const size_t &j, const size_t &i) const {
    return m_g;
}

defs::ham_t HolsteinLadderHam::get_coeff_1101(const size_t &imode, const size_t &j, const size_t &i) const {
    return m_g;
}

defs::ham_t HolsteinLadderHam::get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    auto imode = conn.m_bos.m_cre[0].m_imode;
    auto nocc_frm = onv.m_frm.site_nocc(imode);
    auto occ_fac = std::sqrt(onv.m_bos[imode]+1);//conn.m_bos.occ_fac(onv.m_bos);
    return m_g*nocc_frm*occ_fac;
}

defs::ham_t HolsteinLadderHam::get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    auto imode = conn.m_bos.m_ann[0].m_imode;
    auto nocc_frm = onv.m_frm.site_nocc(imode);
    auto occ_fac = std::sqrt(onv.m_bos[imode]);//conn.m_bos.occ_fac(onv.m_bos);
    return m_g*nocc_frm*occ_fac;
}

defs::ham_t HolsteinLadderHam::get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return 0;
}

defs::ham_t HolsteinLadderHam::get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    return 0;
}