//
// Created by Robert J. Anderson on 12/8/21.
//

#include "GeneralFrmHam.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"
#include "M7_lib/excitgen/frm/Pchb2200.h"

void GeneralFrmHam::log_ints_sym(integrals_1e::syms::Sym sym, bool initial) {
    const str_t context_str = initial ? "initially" : "conflict detected:";
    logging::info("{} assuming {} permutational symmetry for 1e integrals",
              context_str, integrals_1e::syms::name(sym));
}

void GeneralFrmHam::log_ints_sym(integrals_2e::syms::Sym sym, bool initial) {
    const str_t context_str = initial ? "initially" : "conflict detected:";
    logging::info("{} assuming {} permutational symmetry for 2e integrals",
              context_str, integrals_2e::syms::desc(sym));
    const auto equivs = integrals_2e::syms::equivalences(sym);
    if (equivs.size()<2) return;
    logging::info("this storage scheme assumes that {} integrals are equivalent", convert::to_string(equivs));
}


uint_t GeneralFrmHam::make_ints_(IntegralReader *reader, GeneralFrmHam::ints_1e_t *ints_1e, GeneralFrmHam::ints_2e_t *ints_2e) {
    REQUIRE_TRUE(reader, "this method should not be called on a non-reading rank");
    IntegralReader::IterData d;

    while (reader->next(d)){
        if (d.m_ranksig == opsig::c_invalid){
            // non-coefficient entry
            continue;
        }
        if (d.m_ranksig == opsig::c_zero) {
            m_e_core = d.m_value;
            continue;
        }

        auto& rank_contrib = d.m_ranksig == opsig::c_sing ? m_contribs_1100 : m_contribs_2200;
        rank_contrib.set_nonzero(d.m_exsig);

        if (d.m_ranksig == opsig::c_sing) {
            if (!ints_1e->set(d.m_inds[0], d.m_inds[1], d.m_value)) {
                logging::info_("integral < {} | h_core | {} > value {} is at odds with stored value {}",
                               d.m_inds[0], d.m_inds[1], d.m_value, ints_1e->get(d.m_inds[0], d.m_inds[1]));
                return 1ul;
            }
        } else if (d.m_ranksig == opsig::c_doub) {
            // FCIDUMP integral indices are in chemists' ordering
            if (!ints_2e->set(d.m_inds[0], d.m_inds[2], d.m_inds[1], d.m_inds[3], d.m_value)) {
                logging::info_("integral < {} {} | {} {} > value {} is at odds with stored value {}",
                               d.m_inds[0], d.m_inds[2], d.m_inds[1], d.m_inds[3], d.m_value,
                               ints_2e->get(d.m_inds[0], d.m_inds[2], d.m_inds[1], d.m_inds[3]));
                return 2ul;
            }
        } else MPI_ABORT("File reader error");
    }
    return 0ul;
}

GeneralFrmHam::Integrals GeneralFrmHam::make_ints(IntegralReader* reader) {
    if (!m_basis.m_nsite) return {nullptr, nullptr};

    if (reader) {
        m_kramers_attrs.m_conserving_singles = reader->spin_conserving(1);
        m_kramers_attrs.m_conserving_doubles = reader->spin_conserving(2);
    }

    REQUIRE_EQ(m_basis.m_abgrp_map.m_site_irreps.size(),m_basis.m_nsite, "site map size incorrect");

    m_complex_valued = false;
    if (reader) m_complex_valued = reader->complex_valued();
    // get a vector of all reading ranks
    const auto reading_iranks = mpi::filter(bool(reader));
    // use the first one for the definitive value of the complex valued flag
    mpi::bcast(m_complex_valued, reading_iranks.front());

    using namespace ham;
    /*
     * initialize permutational symmetries.
     */
    auto ints_1e = integrals_1e::make<ham_t>(
            m_basis.ncoeff_ind(m_info.m_spin_resolved), integrals_1e::syms::H);
    log_ints_sym(ints_1e->sym(), true);
    REQUIRE_TRUE(ints_1e.get(), "1e integral array object unallocated");
    /*
     * even if the source is complex-valued, it can still have DHR symmetry
     */
    auto ints_2e = integrals_2e::make<ham_t>(m_basis.ncoeff_ind(m_info.m_spin_resolved), m_info.m_init_2e_perm_sym);
    log_ints_sym(ints_2e->sym(), true);
    REQUIRE_TRUE(ints_2e.get(), "2e integral array object unallocated");

    if (reader) m_e_core = reader->ecore();

    uint_t iex_next = ~0ul;
    while (iex_next) {
        if (reader) iex_next = make_ints_(reader, ints_1e.get(), ints_2e.get());
        mpi::bcast(iex_next, reading_iranks.front());
        if (iex_next==1ul) {
            integrals_1e::next_sym_attempt(ints_1e);
            log_ints_sym(ints_1e->sym(), false);
            if (reader) reader->goto_first_1e();
        }
        else if (iex_next==2ul) {
            integrals_2e::next_sym_attempt(ints_2e);
            log_ints_sym(ints_2e->sym(), false);
            if (reader) reader->goto_first_2e();
        }
    }

    mpi::bcast(m_e_core, reading_iranks.front());
    mpi::barrier();
    logging::info("FCIDUMP loading complete.");
    log_data();
    return {std::move(ints_1e), std::move(ints_2e)};
}

GeneralFrmHam::Integrals GeneralFrmHam::make_ints() {
    const str_t fmt = "Reading fermion Hamiltonian coefficients from {} file \"" + m_info.m_fname + "\"...";
    if (m_info.m_impl==FcidumpInfo::CSV) {
        logging::info(fmt, "plain text CSV");

        if (m_info.m_spin_resolved) {
            logging::info("FCIDUMP file is spin resolved (e.g. UHF)");
            logging::info("Following the {} convention for spin-resolved integrals",
                          FcidumpInfo::ur_desc(m_info.m_ur_style));
        }
        else logging::info("FCIDUMP file is not spin resolved (e.g. RHF)");

        if (mpi::on_node_i_am_root()) {
            CsvIntegralReader reader(m_info);
            return make_ints(&reader);
        }
        else return make_ints(nullptr);
    }
    else if (m_info.m_impl==FcidumpInfo::MolcasHDF5) {
        logging::info(fmt, "Molcas HDF5 binary");
        // HDF5 reader needs to be initialized on every rank
        MolcasHdf5IntegralReader reader(m_info);
        return make_ints(&reader);
    }
    return {nullptr, nullptr};
}


GeneralFrmHam::GeneralFrmHam(const FcidumpInfo& info):
        FrmHam({info.m_nsite, {PointGroup(), info.m_orbsym}}), m_info(info), m_ints(make_ints()){
    /*
     * since reading only happens on a subset of ranks, the excitation signatures contributing to each rank were only
     * set on these reading ranks. here they must be synchronized
     */
    m_contribs_1100.bcast();
    m_contribs_2200.bcast();
    /*
     * same for the Kramers symmetry attributes
     */
    m_kramers_attrs.bcast();
}

GeneralFrmHam::GeneralFrmHam(init_opts_t opts):
    GeneralFrmHam(
        FcidumpInfo(
            opts.m_ham.m_fcidump.m_path,
            FcidumpInfo::ur_style(opts.m_ham.m_fcidump.m_unrestrict_style),
            integrals_2e::syms::from_symbol(opts.m_ham.m_fcidump.m_init_2e_perm_sym)
        )
    ){}

ham_t GeneralFrmHam::get_coeff_1100(uint_t a, uint_t i) const {
    if (m_info.m_spin_resolved) return m_ints.m_1e->get(a, i);
    if (m_basis.ispin(i)!=m_basis.ispin(a)) return 0.0;
    return m_ints.m_1e->get(m_basis.isite(a), m_basis.isite(i));
}

ham_t GeneralFrmHam::get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const {
    if (m_info.m_spin_resolved) return m_ints.m_2e->get(a, b, i, j)-m_ints.m_2e->get(a, b, j, i);
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

ham_t GeneralFrmHam::get_element_0000(const field::FrmOnv& onv) const {
    ham_t element = m_e_core;
    auto singles_fn = [&](uint_t i) { element += GeneralFrmHam::get_coeff_1100(i, i); };
    auto doubles_fn = [&](uint_t i, uint_t j) {
        element += GeneralFrmHam::get_coeff_2200(i, j, i, j);
    };
    onv.foreach_setbit_pair(singles_fn, doubles_fn);
    return element;
}

ham_t GeneralFrmHam::get_element_1100(const field::FrmOnv& onv, const conn::FrmOnv& conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), opsig::c_sing, "expected 1100 (aka fermion single) exsig");
    const auto& ann = conn.m_ann[0];
    const auto& cre = conn.m_cre[0];

    ham_t element = GeneralFrmHam::get_coeff_1100(cre, ann);
    auto fn = [&](uint_t ibit) {
        if (ibit != ann) element += GeneralFrmHam::get_coeff_2200(cre, ibit, ann, ibit);
    };
    onv.foreach_setbit(fn);
    return conn.phase(onv) ? -element : element;
}

ham_t GeneralFrmHam::get_element_2200(const field::FrmOnv& onv, const conn::FrmOnv& conn) const {
    DEBUG_ASSERT_EQ(conn.exsig(), opsig::c_doub, "expected 2200 (aka fermion double) exsig");
    const auto element = GeneralFrmHam::get_coeff_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
    return conn.phase(onv) ? -element : element;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(
        PRNG& prng, const conf::Propagator& /*opts*/, const FrmHam& h) {
    excit_gen_list_t list;
    bool any_singles = h.m_contribs_1100.is_nonzero(opsig::c_sing) || h.m_contribs_2200.is_nonzero(opsig::c_sing);
    if (any_singles) list.emplace_front(new exgen::UniformSingles(h, prng));
    bool any_doubles = h.m_contribs_2200.is_nonzero(opsig::c_doub);
    if (any_doubles) list.emplace_front(new exgen::Pchb2200(h, prng));
    return list;
}

HamOpTerm::excit_gen_list_t GeneralFrmHam::make_excit_gens(
        PRNG& prng, const conf::Propagator& opts) const {
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

uint_t GeneralFrmHam::default_nelec() const {
    return m_info.m_nelec;
}

int GeneralFrmHam::default_ms2_value() const {
    return m_info.m_ms2;
}

bool GeneralFrmHam::is_hermitian() const {
    const auto s1 = m_ints.m_1e->sym();
    const auto s2 = m_ints.m_2e->sym();
    const auto h1 = s1==integrals_1e::syms::H;
    const auto h2 = (s2==integrals_2e::syms::DHR) || (s2==integrals_2e::syms::DH) || (s2==integrals_2e::syms::H);
    return h1 && h2;
}
