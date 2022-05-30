//
// Created by Robert J. Anderson on 05/11/2020.
//

#include "FrmBosHam.h"

#include <utility>

FrmBosHam::FrmBosHam(const sys::Basis &basis) :
        m_basis(basis),
        m_contribs_1110(exsig_utils::ex_1110), m_contribs_1101(exsig_utils::ex_1101) {}

void FrmBosHam::log_data() const {
    if (!*this) return;
    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_0010))
        log::info("1110 fermion-coupled boson ladder term has no 0010 contributions");
    if (!m_contribs_1110.is_nonzero(exsig_utils::ex_1110))
        log::info("1110 fermion-coupled boson ladder term has no 1110 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_0001))
        log::info("1101 fermion-coupled boson ladder term has no 0001 contributions");
    if (!m_contribs_1101.is_nonzero(exsig_utils::ex_1101))
        log::info("1101 fermion-coupled boson ladder term has no 1101 contributions");
}