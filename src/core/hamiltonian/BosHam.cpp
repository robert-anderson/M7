//
// Created by rja on 26/07/2021.
//

#include "BosHam.h"
#include "io/BosdumpFileReader.h"


BosHam::BosHam(size_t nmode, size_t nboson):
        m_nmode(nmode), m_nboson(nboson),
        m_contribs_0011(exsig_utils::ex_0011), m_contribs_0022(exsig_utils::ex_0022) {}

size_t BosHam::nci() const {
    return ci_utils::boson_dim(m_nmode, m_nboson, true);
}

void BosHam::log_data() const {
    if (disabled()) return;
    if (!m_contribs_0011.is_nonzero(0ul))
        log::info("1-boson (0011) term has no diagonal (0000) contributions");
    if (!m_contribs_0011.is_nonzero(exsig_utils::ex_0011))
        log::info("1-boson (0011) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(0ul))
        log::info("2-boson (0022) term has no diagonal (0000) contributions");
    if (!m_contribs_0022.is_nonzero(exsig_utils::ex_0011))
        log::info("2-boson (0022) term has no single-excitation (0011) contributions");
    if (!m_contribs_0022.is_nonzero(exsig_utils::ex_0022))
        log::info("2-boson (0022) term has no double-excitation (0022) contributions");
}
