//
// Created by Robert J. Anderson on 27/02/2020.
//

#include "M7_lib/caches/DecodedDeterminants.h"

#include "FrmHam.h"


FrmHam::FrmHam(const sys::frm::Basis& basis):
        m_basis(basis), m_contribs_1100(exsig::ex_single), m_contribs_2200(exsig::ex_double){}

ham_t FrmHam::get_element(const field::FrmOnv &onv) const {
    return get_element_0000(onv);
}

ham_comp_t FrmHam::get_energy(const field::FrmOnv &onv) const {
    auto elem = get_element_0000(onv);
    DEBUG_ASSERT_TRUE(fptol::near_zero(elem), "energies should be purely real");
    return arith::real(elem);
}

ham_t FrmHam::get_element(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    switch (conn.size()) {
        case 0:
            return get_element_0000(onv);
        case 2:
            return get_element_1100(onv, conn);
        case 4:
            return get_element_2200(onv, conn);
        case 6:
            return get_element_3300(onv, conn);
        default:
            return 0.0;
    }
}

void FrmHam::log_data() const {
    if (!*this) return;
    if (!m_contribs_1100.is_nonzero(0ul))
        logging::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(exsig::ex_single))
        logging::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(0ul))
        logging::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(exsig::ex_single))
        logging::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(exsig::ex_double))
        logging::info("2-electron term has no double-excitation contributions");
}