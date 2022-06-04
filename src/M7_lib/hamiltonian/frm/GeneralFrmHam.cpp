//
// Created by Robert J. Anderson on 12/8/21.
//

#include "GeneralFrmHam.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"


GeneralFrmHam::GeneralFrmHam(const FcidumpInfo& info, bool spin_major):
        FrmHam({info.m_nsite, {PointGroup(), info.m_orbsym}, info.m_spin_resolved}),
        m_info(info) {
    if (!m_basis.m_nsite) return;

    REQUIRE_EQ(m_basis.m_abgrp_map.m_site_irreps.size(),m_basis.ncoeff_ind(),"site map size incorrect");

    FcidumpFileReader file_reader(info.m_fname, spin_major);
    m_complex_valued = file_reader.m_complex_valued;

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;
    // assume no contribution cases to be nonzero unless a counterexample is found

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
//            m_int_1.set(inds, value);
        } else if (ranksig == ex_double) {
//            m_int_2.set(inds, value);
        } else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");
    log_data();
}

GeneralFrmHam::GeneralFrmHam(opt_pair_t opts):
        GeneralFrmHam({FortranNamelistReader(opts.m_ham.m_fcidump.m_path)}, opts.m_ham.m_fcidump.m_spin_major) {}

defs::ham_t GeneralFrmHam::get_coeff_1100(size_t a, size_t i) const {
    return m_ints_1e->get(a, i);
}

defs::ham_t GeneralFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    return m_ints_2e->get(a, b, i, j)-m_ints_2e->get(a, b, j, i);
}

defs::ham_t GeneralFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_e_core;
    auto singles_fn = [&](size_t i) { element += m_ints_1e->get(i, i); };
    auto doubles_fn = [&](size_t i, size_t j) {
        element += GeneralFrmHam::get_coeff_2200(i, j, i, j);
    };
    onv.foreach_setbit_pair(singles_fn, doubles_fn);
    return element;
}

defs::ham_t GeneralFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_single, "expected 1100 (aka fermion single) exsig");
    const auto &ann = conn.m_ann[0];
    const auto &cre = conn.m_cre[0];

    defs::ham_t element = m_ints_1e->get(cre, ann);
    auto fn = [&](size_t ibit) {
        if (ibit != ann) element += GeneralFrmHam::get_coeff_2200(cre, ibit, ann, ibit);
    };
    onv.foreach_setbit(fn);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t GeneralFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
    const auto element = GeneralFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(
        PRNG &prng, const conf::Propagator &opts, const FrmHam &h) {
    using namespace exsig_utils;
    excit_gen_list_t list;
    bool any_singles = h.m_contribs_1100.is_nonzero(ex_single) || h.m_contribs_2200.is_nonzero(ex_single);
    if (any_singles) list.emplace_front(new UniformSingles(h, prng));
    bool any_doubles = h.m_contribs_2200.is_nonzero(ex_double);
    if (any_doubles) list.emplace_front(new Pchb2200(h, prng));
    return list;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(
        PRNG &prng, const conf::Propagator &opts) const {
    return make_excit_gens(prng, opts, *this);
}

conn_foreach::base_list_t GeneralFrmHam::make_foreach_iters() const {
    conn_foreach::base_list_t list;
    if (m_kramers_attrs.m_conserving_singles)
        list.emplace_front(new conn_foreach::frm::Ms2Conserve<1>);
    else
        list.emplace_front(new conn_foreach::frm::General<1>);
    if (m_kramers_attrs.m_conserving_doubles)
        list.emplace_front(new conn_foreach::frm::Ms2Conserve<2>);
    else
        list.emplace_front(new conn_foreach::frm::General<2>);
    return list;
}

size_t GeneralFrmHam::default_nelec() const {
    return m_info.m_nelec;
}

int GeneralFrmHam::default_ms2_value() const {
    return m_info.m_ms2;
}
