//
// Created by anderson on 12/8/21.
//

#include "GeneralFrmHam.h"


buffered::FrmOnv GeneralFrmHam::guess_reference(const int &spin_restrict) const {
    buffered::FrmOnv ref({m_nsite, 0ul});
    REQUIRE_EQ(size_t(std::abs(spin_restrict) % 2), m_nelec % 2,
               "Sz quantum number given incompatible with nelec");
    size_t n_spin_0 = (m_nelec + spin_restrict) / 2;
    size_t n_spin_1 = m_nelec - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    DEBUG_ASSERT_EQ(ref.ms2(), spin_restrict, "constructed fermion ONV does not have expected Sz");
    return ref;
}

GeneralFrmHam::GeneralFrmHam(size_t nelec, size_t nsite, bool spin_resolved,
                                         int ms2_restrict, bool complex_valued, defs::inds site_irreps) :
        FrmHam(nelec, nsite, ms2_restrict, complex_valued, site_irreps),
        m_int_1(nsite, spin_resolved), m_int_2(nsite, spin_resolved),
        m_contribs_1100(exsig_utils::ex_single), m_contribs_2200(exsig_utils::ex_double) {
    if (!nsite) return;
    REQUIRE_EQ(m_point_group_map.m_site_irreps.size(), norb_distinct(), "site map size incorrect");
}

GeneralFrmHam::GeneralFrmHam(const FcidumpHeader& header, bool spin_major, int ms2_restrict, int charge) :
        GeneralFrmHam(header.m_nelec - charge, header.m_nsite, header.m_spin_resolved, ms2_restrict,
                           FcidumpFileReader(header.m_fname, spin_major).m_complex_valued, header.m_orbsym) {

    FcidumpFileReader file_reader(header.m_fname, spin_major);

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;

    // assume no contribution cases to be nonzero unless a counterexample is found
    REQUIRE_LE_ALL(m_nelec, 2 * m_nsite, "unphysical number of electrons specified in FCIDUMP file");

    log::info("Reading fermion Hamiltonian coefficients from FCIDUMP file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        if (consts::float_is_zero(value)) continue;
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


defs::ham_t GeneralFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_e_core;
    auto singles_fn = [&](const size_t &i) { element += m_int_1(i, i); };
    auto doubles_fn = [&](const size_t &i, const size_t &j) {
        element += m_int_2.phys_antisym_element(i, j, i, j);
    };
    onv.foreach_pair(singles_fn, doubles_fn);
    return element;
}

defs::ham_t GeneralFrmHam::get_coeff_1100(const size_t &i, const size_t &j) const {
    return m_int_1(i, j);
}

defs::ham_t GeneralFrmHam::get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}
