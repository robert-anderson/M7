//
// Created by anderson on 12/8/21.
//

#include "GeneralFrmHam.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"

buffered::FrmOnv GeneralFrmHam::guess_reference(const int &spin_restrict) const {
    buffered::FrmOnv ref(m_bd.m_nsite);
    REQUIRE_EQ(size_t(std::abs(spin_restrict) % 2), m_nelec % 2,
               "Sz quantum number given incompatible with nelec");
    size_t n_spin_0 = (m_nelec + spin_restrict) / 2;
    size_t n_spin_1 = m_nelec - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    DEBUG_ASSERT_EQ(ref.ms2(), spin_restrict, "constructed fermion ONV does not have expected Sz");
    return ref;
}
/*
 *
    GeneralFrmHam(const FrmBasisData& bd, bool spin_resolved, int ms2_restrict, defs::inds site_irreps = {});

    GeneralFrmHam(const FcidumpHeader& header, bool spin_major, int ms2_restrict, int charge = 0);

    GeneralFrmHam(std::string fname, bool spin_major, int charge = 0):
            GeneralFrmHam(FcidumpHeader(fname), spin_major, charge){}
 */

GeneralFrmHam::GeneralFrmHam(const FrmBasisData &bd, size_t nelec, int ms2_restrict) :
    FrmHam(bd, nelec, ms2_restrict),
    m_int_1(m_bd.m_nsite, m_bd.m_spin_resolved), m_int_2(m_bd.m_nsite, m_bd.m_spin_resolved) {
    if (!m_bd.m_nsite) return;
    REQUIRE_EQ(m_bd.m_abgrp_map.m_site_irreps.size(), norb_distinct(), "site map size incorrect");
}

GeneralFrmHam::GeneralFrmHam(const FcidumpHeader& header, bool spin_major, int ms2_restrict, int charge) :
        GeneralFrmHam({header.m_nsite, {PointGroup(), header.m_orbsym}, header.m_spin_resolved},
                      header.m_nelec - charge, ms2_restrict){

    FcidumpFileReader file_reader(header.m_fname, spin_major);
    m_complex_valued = file_reader.m_complex_valued;

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;

    // assume no contribution cases to be nonzero unless a counterexample is found
    REQUIRE_LE_ALL(m_nelec, m_bd.m_nspinorb,
                   "unphysical number of electrons specified by FCIDUMP file header and charge parameter");

    log::info("Reading fermion Hamiltonian coefficients from FCIDUMP file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::nearly_zero(value)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);

        if (ranksig == 0ul) {
            m_e_core = value;
            continue;
        }

        auto &rank_contrib = ranksig == ex_single ? m_contribs_1100 : m_contribs_2200;
        rank_contrib.set_nonzero(exsig);

        if (ranksig == ex_single) {
            m_int_1.set(inds, value);
        } else if (ranksig == ex_double) {
            m_int_2.set(inds, value);
        } else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");
    log_data();
}


defs::ham_t GeneralFrmHam::get_coeff_1100(size_t a, size_t i) const {
    return m_int_1(a, i);
}

defs::ham_t GeneralFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    return m_int_2.phys_antisym_element(a, b, i, j);
}

defs::ham_t GeneralFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_e_core;
    auto singles_fn = [&](size_t i) { element += m_int_1(i, i); };
    auto doubles_fn = [&](size_t i, size_t j) {
        element += m_int_2.phys_antisym_element(i, j, i, j);
    };
    onv.foreach_setbit_pair(singles_fn, doubles_fn);
    return element;
}

defs::ham_t GeneralFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_single, "expected 1100 (aka fermion single) exsig");
    const auto &ann = conn.m_ann[0];
    const auto &cre = conn.m_cre[0];

    defs::ham_t element = m_int_1(cre, ann);
    auto fn = [&](size_t ibit) {
        if (ibit != ann) element += m_int_2.phys_antisym_element(cre, ibit, ann, ibit);
    };
    onv.foreach_setbit(fn);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t GeneralFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
    const auto element = m_int_2.phys_antisym_element(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(PRNG &prng, const fciqmc_config::Propagator& opts) const {
    using namespace exsig_utils;
    excit_gen_list_t list;
    bool any_singles = m_contribs_1100.is_nonzero(ex_single) || m_contribs_2200.is_nonzero(ex_single);
    if (any_singles) list.emplace_front(new UniformSingles(*this, prng));
    bool any_doubles = m_contribs_2200.is_nonzero(ex_double);
    if (any_doubles) list.emplace_front(new Pchb2200(*this, prng));
    return list;
}

conn_foreach::base_list_t GeneralFrmHam::make_foreach_iters() const {
    conn_foreach::base_list_t list;
    if (m_kramers_attrs.m_conserving_singles)
        list.emplace_front(new conn_foreach::frm::Ms2Conserve<1>(m_bd));
    else
        list.emplace_front(new conn_foreach::frm::General<1>(m_bd));
    if (m_kramers_attrs.m_conserving_doubles)
        list.emplace_front(new conn_foreach::frm::Ms2Conserve<2>(m_bd));
    else
        list.emplace_front(new conn_foreach::frm::General<2>(m_bd));
    return list;
}