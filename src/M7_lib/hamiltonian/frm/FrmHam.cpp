//
// Created by rja on 27/02/2020.
//

#include "M7_lib/caches/DecodedDeterminants.h"

#include "FrmHam.h"

FrmHam::FrmHam(const sys::frm::Basis& basis):
        m_basis(basis), m_contribs_1100(exsig_utils::ex_single), m_contribs_2200(exsig_utils::ex_double) {}

defs::ham_t FrmHam::get_element(const field::FrmOnv &onv) const {
    return get_element_0000(onv);
}

defs::ham_comp_t FrmHam::get_energy(const field::FrmOnv &onv) const {
    return consts::real(get_element_0000(onv));
}

defs::ham_t FrmHam::get_element(const field::FrmOnv &ket, const conn::FrmOnv &conn) const {
    switch (conn.size()) {
        case 0:
            return get_element_0000(ket);
        case 2:
            return get_element_1100(ket, conn);
        case 4:
            return get_element_2200(ket, conn);
        default:
            return 0.0;
    }
}

void FrmHam::log_data() const {
    if (disabled()) return;
    if (!m_contribs_1100.is_nonzero(0ul))
        log::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(exsig_utils::ex_single))
        log::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(0ul))
        log::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_single))
        log::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_double))
        log::info("2-electron term has no double-excitation contributions");
}
