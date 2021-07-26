//
// Created by rja on 24/05/2021.
//

#include <src/defs.h>
#include "SymmetryHelpers.h"
#include "Hamiltonian.h"

#if 0
ham_sym_helpers::Fermion::Fermion(const FermionHamiltonian &ham) :
        m_ham(ham), m_nsite(ham.nsite()), m_conn_work(m_nsite),
        m_occ_work(m_nsite), m_vac_work(m_nsite), m_onv_work(m_nsite) {}

defs::ham_t ham_sym_helpers::Fermion::get_element(const conn::Antisym<0> &conn) const {
    return m_ham.get_element(conn);
}

defs::ham_comp_t ham_sym_helpers::Fermion::get_energy(const fields::Onv<0> &onv) const {
    m_conn_work.connect(onv, onv);
    return consts::real(get_element(m_conn_work));
}

bool ham_sym_helpers::Fermion::update_helement(bool get_h, bool h_nonzero_only) const {
    if (get_h || h_nonzero_only) m_helement_work = get_element(m_conn_work);
    else m_helement_work = 0.0;
    if (h_nonzero_only) return !consts::float_is_zero(m_helement_work);
    return true;
}

void ham_sym_helpers::Fermion::perform_diagonal(const fields::Onv<0> &src_onv,
                                                const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
                                                bool h_nonzero_only) const {
    m_conn_work.zero();
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);
}

void ham_sym_helpers::Fermion::perform_single(const fields::Onv<0> &src_onv, const size_t &occ, const size_t &vac,
                                              const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
                                              bool h_nonzero_only) const {
    m_conn_work.zero();
    m_conn_work.add(occ, vac);
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);
}

void ham_sym_helpers::Fermion::perform_double(const fields::Onv<0> &src_onv, const size_t &occ1, const size_t &occ2,
                                              const size_t &vac1, const size_t &vac2,
                                              const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
                                              bool h_nonzero_only) const {
    m_conn_work.zero();
    m_conn_work.add(occ1, occ2, vac1, vac2);
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);
}

void ham_sym_helpers::Fermion::foreach_connection_singles(
        const fields::Onv<0> &src_onv, const defs::inds &occs,
        const defs::inds &vacs,
        const ham_sym_helpers::Fermion::body_fn_t &body,
        bool get_h, bool h_nonzero_only) const {
    for (size_t iocc = 0ul; iocc < occs.size(); ++iocc) {
        auto &occ = occs[iocc];
        for (size_t ivac = 0ul; ivac < vacs.size(); ++ivac) {
            auto &vac = vacs[ivac];
            perform_single(src_onv, occ, vac, body, get_h, h_nonzero_only);
        }
    }
}

void ham_sym_helpers::Fermion::foreach_connection_subset_doubles(
        const fields::Onv<0> &src_onv, const defs::inds &occs1,
        const defs::inds &occs2, const defs::inds &vacs1,
        const defs::inds &vacs2,
        const ham_sym_helpers::Fermion::body_fn_t &body,
        bool get_h,
        bool h_nonzero_only) const {
    for (size_t iocc1 = 0ul; iocc1 < occs1.size(); ++iocc1) {
        auto &occ1 = occs1[iocc1];
        for (size_t ivac1 = 0ul; ivac1 < vacs1.size(); ++ivac1) {
            auto &vac1 = vacs1[ivac1];
            size_t iocc2_start = (&occs1 == &occs2) ? iocc1 + 1 : 0ul;
            for (size_t iocc2 = iocc2_start; iocc2 < occs2.size(); ++iocc2) {
                auto &occ2 = occs2[iocc2];
                size_t ivac2_start = (&vacs1 == &vacs2) ? ivac1 + 1 : 0ul;
                for (size_t ivac2 = ivac2_start; ivac2 < vacs2.size(); ++ivac2) {
                    auto &vac2 = vacs2[ivac2];
                    perform_double(src_onv, occ1, occ2, vac1, vac2, body, get_h, h_nonzero_only);
                }
            }
        }
    }
}

void ham_sym_helpers::Fermion::foreach_connection_subset_doubles(
        const fields::Onv<0> &src_onv, const defs::inds &occs,
        const defs::inds &vacs,
        const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
        bool h_nonzero_only) const {
    foreach_connection_subset_doubles(src_onv, occs, vacs, occs, vacs, body, get_h, h_nonzero_only);
}

void ham_sym_helpers::Fermion::foreach_connection_subset(
        const fields::Onv<0> &src_onv, const defs::inds &occs1,
        const defs::inds &occs2, const defs::inds &vacs1,
        const defs::inds &vacs2,
        const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
        bool h_nonzero_only) const {
    for (size_t iocc1 = 0ul; iocc1 < occs1.size(); ++iocc1) {
        auto &occ1 = occs1[iocc1];
        for (size_t ivac1 = 0ul; ivac1 < vacs1.size(); ++ivac1) {
            auto &vac1 = vacs1[ivac1];
            perform_single(src_onv, occ1, vac1, body, get_h, h_nonzero_only);
            size_t iocc2_start = (&occs1 == &occs2) ? iocc1 + 1 : 0ul;
            for (size_t iocc2 = iocc2_start; iocc2 < occs2.size(); ++iocc2) {
                auto &occ2 = occs2[iocc2];
                size_t ivac2_start = (&vacs1 == &vacs2) ? ivac1 + 1 : 0ul;
                for (size_t ivac2 = ivac2_start; ivac2 < vacs2.size(); ++ivac2) {
                    auto &vac2 = vacs2[ivac2];
                    perform_double(src_onv, occ1, occ2, vac1, vac2, body, get_h, h_nonzero_only);
                }
            }
        }
    }
}

void ham_sym_helpers::Fermion::foreach_connection_subset(
        const fields::Onv<0> &src_onv, const defs::inds &occs,
        const defs::inds &vacs,
        const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
        bool h_nonzero_only) const {
    foreach_connection_subset(src_onv, occs, occs, vacs, vacs, body, get_h, h_nonzero_only);
}

void ham_sym_helpers::Fermion::foreach_connection(
        const fields::Onv<0> &src_onv,
        const ham_sym_helpers::Fermion::body_fn_t &body, bool get_h,
        bool h_nonzero_only, bool include_diagonal) const {
    m_helement_work = 0.0;
    m_occ_work.update(src_onv);
    m_vac_work.update(src_onv);
    if (include_diagonal) perform_diagonal(src_onv, body, get_h, h_nonzero_only);
    foreach_connection_subset(src_onv, m_occ_work.inds(), m_vac_work.inds(), body, get_h, h_nonzero_only);
}


ham_sym_helpers::FermiBos::FermiBos(const FermiBosHamiltonian &ham) :
        m_ham(ham), m_nsite(ham.nsite()), m_conn_work(m_nsite),
        m_occ_work(m_nsite), m_vac_work(m_nsite), m_onv_work(m_nsite) {}



defs::ham_t ham_sym_helpers::FermiBos::get_element(const conn::Antisym<1> &conn) const {
    return m_ham.get_element(conn);
}

defs::ham_comp_t ham_sym_helpers::FermiBos::get_energy(const fields::Onv<1> &onv) const {
    m_conn_work.connect(onv, onv);
    return consts::real(get_element(m_conn_work));
}


void ham_sym_helpers::FermiBos::foreach_connection(
        const fields::Onv<1> &src_onv,
        const ham_sym_helpers::FermiBos::body_fn_t &body, bool get_h,
        bool h_nonzero_only, bool include_diagonal) const {

    m_conn_work.zero();
    m_conn_work.apply(src_onv, m_onv_work);
    if (update_helement(get_h, h_nonzero_only)) body(m_conn_work, m_onv_work, m_helement_work);

    m_conn_work.m_bonvconn.zero();
    m_onv_work.m_bos = src_onv.m_bos;

    //static_cast<const FermionHamiltonian&>(m_ham).foreach_connection()
//    auto frm_body = [&](const conn::Antisym<0> &, const fields::Onv<0> &, const defs::ham_t &){
//        body(m_conn_work, m_onv_work, m_helement_work);
//    };
    // do all purely fermionic connections:
    //m_ham.foreach_connection(src_onv.m_frm, frm_body, get_h, h_nonzero_only, include_diagonal);
}
#endif