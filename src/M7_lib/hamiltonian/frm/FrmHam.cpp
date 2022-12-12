//
// Created by Robert J. Anderson on 27/02/2020.
//

#include "M7_lib/caches/DecodedDeterminants.h"

#include "FrmHam.h"


FrmHam::FrmHam(const sys::frm::Basis& basis):
        m_basis(basis), m_contribs_1100(opsig::c_sing), m_contribs_2200(opsig::c_doub){}

ham_t FrmHam::get_element(const field::FrmOnv &onv) const {
    return get_element_0000(onv);
}

ham_comp_t FrmHam::get_energy(const field::FrmOnv &onv) const {
    auto elem = get_element_0000(onv);
    DEBUG_ASSERT_TRUE(fptol::near_real(elem), "energies should be purely real");
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
    if (!m_contribs_1100.is_nonzero(opsig::c_zero))
        logging::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(opsig::c_sing))
        logging::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(opsig::c_zero))
        logging::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(opsig::c_sing))
        logging::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(opsig::c_doub))
        logging::info("2-electron term has no double-excitation contributions");
}

sys::frm::Electrons FrmHam::electrons(const conf::Particles& p) const {
    uint_t nelec = p.m_nelec.m_value;
    int ms2_value = p.m_ms2.m_value;
    // if either of these are not set, fallback to defaults.
    if (!nelec) nelec = default_nelec();
    if (ms2_value==sys::frm::c_undefined_ms2) ms2_value = default_ms2_value();
    if (ms2_value==sys::frm::c_undefined_ms2) ms2_value = sys::frm::Ms2::lowest_value(nelec);
    return {nelec, ms2_value};
}
