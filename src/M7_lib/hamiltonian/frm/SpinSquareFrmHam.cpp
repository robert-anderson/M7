//
// Created by Oskar Weser on 17/3/22.
//

#include "SpinSquareFrmHam.h"

SpinSquareFrmHam::SpinSquareFrmHam(const sys::frm::Basis &hs):
    FrmHam(hs), m_sz_term(0.25 * m_hs.m_ms2 * (m_hs.m_ms2 - 2)){
    REQUIRE_NE(m_hs.m_ms2, ~0, "spin square operator requires a Hilbert space restricted to a given 2*Ms sector");
}

SpinSquareFrmHam::SpinSquareFrmHam(const FrmHam &h) : SpinSquareFrmHam(h.m_hs){
    REQUIRE_TRUE(h.m_kramers_attrs.conserving(),
                 "spin square operator inconsistent with Kramers non-conservation");
}

defs::ham_t SpinSquareFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    uint n_os_a = 0;
    auto count_n_os_a = [&](size_t i) {
        if (i < m_hs.m_sites) {
            n_os_a += onv.get({1, i});
        };
    };
    onv.foreach_setbit(count_n_os_a);
    return m_sz_term + n_os_a;
}

defs::ham_t SpinSquareFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    return 0.0;
}


defs::ham_t SpinSquareFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    const auto element = this->get_coeff_2200(
                conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t SpinSquareFrmHam::get_coeff_1100(size_t a, size_t i) const {
    return 0.0;
}

defs::ham_t SpinSquareFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    const auto& sites = m_hs.m_sites;
    // We have to determine if it is an exchange.
    return sites.ispin(a) != sites.ispin(b)
        && sites.ispin(i) != sites.ispin(j)
        && sites.isite(a) != sites.isite(b)
        && sites.isite(i) != sites.isite(j)
        && sites.isite(a) == sites.isite(j)
        && sites.isite(b) == sites.isite(i);
}