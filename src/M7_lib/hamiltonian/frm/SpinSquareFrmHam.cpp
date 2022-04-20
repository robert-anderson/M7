//
// Created by Oskar Weser on 17/3/22.
//

#include "SpinSquareFrmHam.h"


SpinSquareFrmHam::SpinSquareFrmHam(size_t nsite, const FrmHilbertData& hd): FrmHam({nsite}, hd){}

SpinSquareFrmHam::SpinSquareFrmHam(const FrmHam &in_ham) : FrmHam(in_ham){
    REQUIRE_TRUE(in_ham.m_kramers_attrs.conserving(),
                 "spin square operator inconsistent with Kramers non-conservation");
}

defs::ham_t SpinSquareFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    uint n_os_a = 0;
    auto count_n_OS_a = [&](size_t i) {
        if (i < m_bd.m_nsite) {
            n_os_a += onv.get({1, i});
        };
    };
    onv.foreach_setbit(count_n_OS_a);
    return m_Sz_term + n_os_a;
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
    auto get_spin = [this](size_t k)
    {
        return field::FrmOnv::ispin(k, m_bd.m_nsite);
    };
    auto spatial_orbital = [this](size_t k)
    {
        return field::FrmOnv::isite(k, m_bd.m_nsite);
    };
    // We have to determine if it is an exchange.
    return get_spin(a) != get_spin(b)
        && get_spin(i) != get_spin(j)
        && spatial_orbital(a) != spatial_orbital(b)
        && spatial_orbital(i) != spatial_orbital(j)
        && spatial_orbital(a) == spatial_orbital(j)
        && spatial_orbital(b) == spatial_orbital(i);
}
