//
// Created by Robert J. Anderson on 05/11/2020.
//

#include "FrmBosHam.h"

#include <utility>

FrmBosHam::FrmBosHam(sys::Basis basis):
        m_basis(std::move(basis)), m_contribs_1110(opsig::c_1110), m_contribs_1101(opsig::c_1101) {}

void FrmBosHam::log_data() const {
    if (!*this) return;
    if (!m_contribs_1110.is_nonzero(opsig::c_0010))
        logging::info("1110 fermion-coupled boson ladder term has no 0010 contributions");
    if (!m_contribs_1110.is_nonzero(opsig::c_1110))
        logging::info("1110 fermion-coupled boson ladder term has no 1110 contributions");
    if (!m_contribs_1101.is_nonzero(opsig::c_0001))
        logging::info("1101 fermion-coupled boson ladder term has no 0001 contributions");
    if (!m_contribs_1101.is_nonzero(opsig::c_1101))
        logging::info("1101 fermion-coupled boson ladder term has no 1101 contributions");
}