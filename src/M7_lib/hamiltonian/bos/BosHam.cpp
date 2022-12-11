//
// Created by Robert J. Anderson on 26/07/2021.
//

#include "M7_lib/io/BosdumpFileReader.h"

#include "BosHam.h"


void BosHam::log_data() const {
    if (*this) return;
    if (!m_contribs_0011.is_nonzero(opsig::c_0000))
        logging::info("1-boson (0011) term has no diagonal (0000) contributions");
    if (!m_contribs_0011.is_nonzero(opsig::c_0011))
        logging::info("1-boson (0011) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(opsig::c_0000))
        logging::info("2-boson (0022) term has no diagonal (0000) contributions");
    if (!m_contribs_0022.is_nonzero(opsig::c_0011))
        logging::info("2-boson (0022) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(opsig::c_0022))
        logging::info("2-boson (0022) term has no double-excitation (0022) contributions");
}
