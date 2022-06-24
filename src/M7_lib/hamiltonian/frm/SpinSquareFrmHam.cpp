//
// Created by Oskar Weser on 17/3/22.
//

#include "SpinSquareFrmHam.h"

SpinSquareFrmHam::SpinSquareFrmHam(const sys::frm::Basis& basis, const sys::frm::Electrons& elecs):
        FrmHam(basis), ElecSpecTerm(elecs), m_sz_term(0.25 * elecs.m_ms2 * (elecs.m_ms2 - 2)){
    REQUIRE_TRUE(elecs.m_ms2.conserve(),
                 "spin square operator requires a Hilbert space restricted to a given 2*Ms sector");
}

defs::ham_t SpinSquareFrmHam::get_element_0000(const field::FrmOnv& onv) const {
    return m_sz_term + onv.nopen_shell_beta();
}

defs::ham_t SpinSquareFrmHam::get_element_2200(const field::FrmOnv& onv, const conn::FrmOnv& conn) const {
    const auto element = this->get_coeff_2200(
                conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t SpinSquareFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    // We have to determine if it is an exchange.
    return m_basis.ispin(a) != m_basis.ispin(b)
        && m_basis.ispin(i) != m_basis.ispin(j)
        && m_basis.isite(a) != m_basis.isite(b)
        && m_basis.isite(i) != m_basis.isite(j)
        && m_basis.isite(a) == m_basis.isite(j)
        && m_basis.isite(b) == m_basis.isite(i);
}