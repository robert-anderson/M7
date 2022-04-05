//
// Created by rja on 05/11/2020.
//

#include "FrmBosHam.h"

FrmBosHam::FrmBosHam(const BasisData &bd, size_t nboson_max) :
        m_bd(bd), m_nboson_max(nboson_max),
        m_contribs_0010(exsig_utils::ex_0010), m_contribs_0001(exsig_utils::ex_0001),
        m_contribs_1110(exsig_utils::ex_1110), m_contribs_1101(exsig_utils::ex_1101) {
    if (!m_nboson_max || !(m_bd.m_nsite || m_bd.m_nmode)) return;
    REQUIRE_EQ(m_bd.m_nsite == 0, m_bd.m_nmode == 0,
               "if the number of sites is non-zero, so also must be the number of boson modes. ");
    log_data();
}

void FrmBosHam::log_data() const {
    if (disabled()) return;
    if (!m_contribs_0010.is_nonzero(exsig_utils::ex_0010))
        log::info("0010 uncoupled boson ladder hamiltonian term has no contributions");
    if (!m_contribs_0001.is_nonzero(exsig_utils::ex_0001))
        log::info("0001 uncoupled boson ladder hamiltonian term has no contributions");

    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_0010))
        log::info("1110 fermion-coupled boson ladder term has no 0010 contributions");
    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_1110))
        log::info("1110 fermion-coupled boson ladder term has no 1110 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_0001))
        log::info("1101 fermion-coupled boson ladder term has no 0001 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_1101))
        log::info("1101 fermion-coupled boson ladder term has no 1101 contributions");
}