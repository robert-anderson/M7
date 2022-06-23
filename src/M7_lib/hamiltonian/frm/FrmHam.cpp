//
// Created by Robert J. Anderson on 27/02/2020.
//

#include "M7_lib/caches/DecodedDeterminants.h"

#include "FrmHam.h"


FrmHam::FrmHam(const sys::frm::Basis& basis):
        m_basis(basis), m_contribs_1100(utils::exsig::ex_single), m_contribs_2200(utils::exsig::ex_double),
        m_work_conn({99, 99}){}

defs::ham_t FrmHam::get_element(const field::FrmOnv &onv) const {
    return get_element_0000(onv);
}

defs::ham_comp_t FrmHam::get_energy(const field::FrmOnv &onv) const {
    auto elem = get_element_0000(onv);
    DEBUG_ASSERT_TRUE(fptol::numeric_real(elem), "energies should be purely real");
    return arith::real(elem);
}

defs::ham_t FrmHam::get_element(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
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
        log::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(utils::exsig::ex_single))
        log::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(0ul))
        log::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(utils::exsig::ex_single))
        log::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(utils::exsig::ex_double))
        log::info("2-electron term has no double-excitation contributions");
}