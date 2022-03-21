//
// Created by Oskar Weser on 17/3/22.
//

#include "SpinHam.h"


defs::ham_t SpinHam::get_element_0000(const field::FrmOnv &onv) const {
    // TODO(@Oskar): make it a const member
    const auto Sz_term = 0.25 * m_ms2_restrict * (m_ms2_restrict - 2);
    uint n_os_a = 0;
    auto count_n_OS_a = [&](size_t i) {
        if (i < m_nsite) {
            n_os_a += onv.get(i + m_nsite);
        };
    };
    onv.foreach(count_n_OS_a);
    return Sz_term + n_os_a;
}


defs::ham_t SpinHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    return 0.0;
}


defs::ham_t SpinHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    const auto element = this->get_coeff_2200(
                conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t SpinHam::get_coeff_1100(size_t i, size_t j) const {
    return 0.0;
}

defs::ham_t SpinHam::get_coeff_2200(size_t i, size_t j, size_t k, size_t l) const {
    auto get_spin = [this](size_t i)
    {
        return field::FrmOnv::ispin(i, m_nsite);
    };
    auto spatial_orbital = [this](size_t i)
    {
        return field::FrmOnv::isite(i, m_nsite);
    };
    return (get_spin(i) == get_spin(j)
            || get_spin(k) == get_spin(l)
            || spatial_orbital(i) == spatial_orbital(j)
            || spatial_orbital(k) == spatial_orbital(l)
            || (spatial_orbital(i) == spatial_orbital(k)
             && spatial_orbital(j) != spatial_orbital(l))
            || (spatial_orbital(i) == spatial_orbital(l)
             && spatial_orbital(j) != spatial_orbital(k))) ? 0.0 : 1.0;
}
