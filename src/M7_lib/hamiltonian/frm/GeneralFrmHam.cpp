//
// Created by Robert J. Anderson on 12/8/21.
//

#include "GeneralFrmHam.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"

void GeneralFrmHam::log_ints_sym(integrals_1e::syms::Sym sym, bool initial) {
    const std::string context_str = initial ? "initially" : "conflict detected:";
    log::info("{} assuming {} permutational symmetry for 1e integrals",
              context_str, integrals_1e::syms::name(sym));
}

void GeneralFrmHam::log_ints_sym(integrals_2e::syms::Sym sym, bool initial) {
    const std::string context_str = initial ? "initially" : "conflict detected:";
    log::info("{} assuming {} permutational symmetry for 2e integrals",
              context_str, integrals_2e::syms::name(sym));
    const auto equivs = integrals_2e::syms::equivalences(sym);
    if (equivs.size()<2) return;
    log::info("this storage scheme assumes that {} integrals are equivalent", utils::convert::to_string(equivs));
}

GeneralFrmHam::Integrals GeneralFrmHam::make_ints(const FcidumpInfo &info, bool spin_major) {
    if (!m_basis.m_nsite) return {nullptr, nullptr};

    REQUIRE_EQ(m_basis.m_abgrp_map.m_site_irreps.size(),m_basis.ncoeff_ind(),"site map size incorrect");

    FcidumpFileReader file_reader(info.m_fname, spin_major);
    m_complex_valued = file_reader.m_complex_valued;

    using namespace ham_data;
    defs::inds_t inds(4);
    defs::ham_t value;

    /*
     * initialize permutational symmetries.
     * if the source is complex-valued, then it cannot have DHR symmetry (and still represent a physical Hamiltonian)
     */
    auto ints_1e = integrals_1e::make<defs::ham_t>(m_basis.ncoeff_ind(), integrals_1e::syms::H);
    log_ints_sym(ints_1e->sym(), true);
    REQUIRE_TRUE(ints_1e.get(), "1e integral array object unallocated");
    auto ints_2e = integrals_2e::make<defs::ham_t>(m_basis.ncoeff_ind(), m_complex_valued ? integrals_2e::syms::DR : integrals_2e::syms::DHR);
    log_ints_sym(ints_2e->sym(), true);
    REQUIRE_TRUE(ints_2e.get(), "2e integral array object unallocated");

    size_t iline_first_2e = ~0ul;
    size_t iline_first_1e = ~0ul;

    log::info("Reading fermion Hamiltonian coefficients from FCIDUMP file \"" + file_reader.m_fname + "\"...");
    size_t iline = ~0ul;
    while (file_reader.next(inds, value)) {
        ++iline;
        if (consts::nearly_zero(value, defs::helem_tol)) continue;
        auto ranksig = file_reader.ranksig(inds);
        auto exsig = file_reader.exsig(inds, ranksig);

        if (ranksig == 0ul) {
            m_e_core = value;
            continue;
        }

        auto &rank_contrib = ranksig == ex_single ? m_contribs_1100 : m_contribs_2200;
        rank_contrib.set_nonzero(exsig);

        if (ranksig == ex_single) {
            if (iline_first_1e==~0ul) iline_first_1e = iline;
            bool success = false;
            while (!success) {
                success = ints_1e->set(inds[0], inds[1], value);
                if (!success) {
                    integrals_1e::next_sym_attempt(ints_1e);
                    file_reader.reset(iline_first_1e);
                    iline = iline_first_1e-1;
                    log_ints_sym(ints_1e->sym(), false);
                }
            }
        } else if (ranksig == ex_double) {
            if (iline_first_2e==~0ul) iline_first_2e = iline;
            bool success = false;
            while (!success) {
                // FCIDUMP integral indices are in chemists' ordering
                success = ints_2e->set(inds[0], inds[2], inds[1], inds[3], value);
                if (!success) {
                    integrals_2e::next_sym_attempt(ints_2e);
                    file_reader.reset(iline_first_2e);
                    iline = iline_first_2e-1;
                    log_ints_sym(ints_2e->sym(), false);
                }
            }

        } else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");
    log_data();
    return {std::move(ints_1e), std::move(ints_2e)};
}



GeneralFrmHam::GeneralFrmHam(const FcidumpInfo& info, bool spin_major):
        FrmHam({info.m_nsite, {PointGroup(), info.m_orbsym}, info.m_spin_resolved}),
        m_info(info), m_ints(make_ints(m_info, spin_major)){
}

GeneralFrmHam::GeneralFrmHam(opt_pair_t opts):
        GeneralFrmHam({FortranNamelistReader(opts.m_ham.m_fcidump.m_path)},
                      opts.m_ham.m_fcidump.m_spin_major) {}

defs::ham_t GeneralFrmHam::get_coeff_1100(size_t a, size_t i) const {
    if (m_basis.m_spin_resolved) return m_ints.m_1e->get(a, i);
    if (m_basis.ispin(i)!=m_basis.ispin(a)) return 0.0;
    return m_ints.m_1e->get(m_basis.isite(a), m_basis.isite(i));
}

defs::ham_t GeneralFrmHam::get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {
    if (m_basis.m_spin_resolved) return m_ints.m_2e->get(a, b, i, j)-m_ints.m_2e->get(a, b, j, i);
    const auto aspin = m_basis.ispin(a);
    const auto bspin = m_basis.ispin(b);
    const auto ispin = m_basis.ispin(i);
    const auto jspin = m_basis.ispin(j);
    if ((aspin + bspin) != (ispin + jspin)) return 0.0;
    a = m_basis.isite(a);
    b = m_basis.isite(b);
    i = m_basis.isite(i);
    j = m_basis.isite(j);
    if (ispin==jspin) return m_ints.m_2e->get(a, b, i, j)-m_ints.m_2e->get(a, b, j, i);
    if (ispin == aspin) return m_ints.m_2e->get(a, b, i, j);
    return -m_ints.m_2e->get(a, b, j, i);
}

defs::ham_t GeneralFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_e_core;
    auto singles_fn = [&](size_t i) { element += GeneralFrmHam::get_coeff_1100(i, i); };
    auto doubles_fn = [&](size_t i, size_t j) {
        element += GeneralFrmHam::get_coeff_2200(i, j, i, j);
    };
    onv.foreach_setbit_pair(singles_fn, doubles_fn);
    return element;
}

defs::ham_t GeneralFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), utils::exsig::ex_single, "expected 1100 (aka fermion single) exsig");
    const auto &ann = conn.m_ann[0];
    const auto &cre = conn.m_cre[0];

    defs::ham_t element = GeneralFrmHam::get_coeff_1100(cre, ann);
    auto fn = [&](size_t ibit) {
        if (ibit != ann) element += GeneralFrmHam::get_coeff_2200(cre, ibit, ann, ibit);
    };
    onv.foreach_setbit(fn);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t GeneralFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), utils::exsig::ex_double, "expected 2200 (aka fermion double) exsig");
    const auto element = GeneralFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(
        PRNG &prng, const conf::Propagator &opts, const FrmHam &h) {
    using namespace utils::exsig;
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