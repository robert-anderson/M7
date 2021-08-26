//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "src/core/basis/DecodedDeterminants.h"
#include "src/core/parallel/RankAllocator.h"
#include "FermionHamiltonian.h"


buffered::FrmOnv FermionHamiltonian::guess_reference(const int &spin_restrict) const {
    buffered::FrmOnv ref(m_nsite);
    REQUIRE_EQ(size_t(std::abs(spin_restrict) % 2), m_nelec % 2,
               "Sz quantum number given incompatible with nelec");
    size_t n_spin_0 = (m_nelec + spin_restrict) / 2;
    size_t n_spin_1 = m_nelec - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    DEBUG_ASSERT_EQ(ref.spin(), spin_restrict, "constructed fermion ONV does not have expected Sz");
    return ref;
}

FermionHamiltonian::FermionHamiltonian(size_t nelec, size_t nsite, bool complex_valued,
                                       bool spin_resolved, defs::inds site_irreps) :
        m_nelec(nelec), m_nsite(nsite), m_complex_valued(complex_valued),
        m_point_group_map(PointGroup(), site_irreps.empty() ? defs::inds(nsite, 0ul) : site_irreps),
        m_int_1(nsite, spin_resolved), m_int_2(nsite, spin_resolved),
        m_contribs_1100(exsig_utils::ex_single), m_contribs_2200(exsig_utils::ex_double){

    DEBUG_ASSERT_EQ(m_point_group_map.m_site_irreps.size(), nsite, "site map size incorrect");
}


FermionHamiltonian::FermionHamiltonian(const FcidumpFileReader &file_reader) :
        FermionHamiltonian(file_reader.m_nelec, file_reader.m_nspatorb,
                           file_reader.m_complex_valued, file_reader.m_spin_resolved, file_reader.m_orbsym) {

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

        if (ranksig==0ul) {
            m_int_0 = value;
            continue;
        }

        auto& rank_contrib = ranksig==ex_single ? m_contribs_1100 : m_contribs_2200;
        rank_contrib.set_nonzero(exsig);

        if (ranksig==ex_single) {
            m_int_1.set(inds, value);
            m_model_attrs.nonzero(m_nsite, inds[0], inds[1]);
        }
        else if (ranksig==ex_double) {
            m_int_2.set(inds, value);
            m_model_attrs.nonzero(m_nsite, inds[0], inds[1], inds[2], inds[3]);
        }
        else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");
    log_data();
}

FermionHamiltonian::FermionHamiltonian(std::string fname, bool spin_major) :
        FermionHamiltonian(FcidumpFileReader(fname, spin_major)) {}


defs::ham_t FermionHamiltonian::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_int_0;
    auto singles_fn = [&](const size_t &i) { element += m_int_1(i, i); };
    auto doubles_fn = [&](const size_t &i, const size_t &j) { element += m_int_2.phys_antisym_element(i, j, i, j); };
    onv.foreach_pair(singles_fn, doubles_fn);
    return element;
}


defs::ham_t
FermionHamiltonian::get_element_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}

void FermionHamiltonian::set_hf_mbf(field::FrmOnv &onv, int spin) const {
    auto nalpha = ci_utils::nalpha(m_nelec, spin);
    auto nbeta = ci_utils::nbeta(m_nelec, spin);
    DEBUG_ASSERT_EQ(nalpha + nbeta, m_nelec, "inconsistent na, nb, nelec");
    onv.zero();
    for (size_t i = 0ul; i < nalpha; ++i) onv.set({0, i});
    for (size_t i = 0ul; i < nbeta; ++i) onv.set({1, i});
}

void FermionHamiltonian::set_afm_mbf(field::FrmOnv &onv, bool alpha_first) const {
    onv.zero();
    size_t ispin = !alpha_first;
    for (size_t isite=0ul; isite<m_nsite; ++isite){
        onv.set({ispin, isite});
        ispin = !ispin;
    }
}

void FermionHamiltonian::log_data() const {
    if (!m_contribs_1100.is_nonzero(0ul))
        log::info("1-electron term has no diagonal contributions");
    if (!m_contribs_1100.is_nonzero(exsig_utils::ex_single))
        log::info("1-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(0ul))
        log::info("2-electron term has no diagonal contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_single))
        log::info("2-electron term has no single-excitation contributions");
    if (!m_contribs_2200.is_nonzero(exsig_utils::ex_double))
        log::info("2-electron term has no double-excitation contributions");
    if (m_model_attrs.m_nn_only_singles)
        log::info("single-excitation contributions to 1-electron term are nearest-neighbor only");
    else if (m_model_attrs.m_nnp_only_singles)
        log::info("single-excitation contributions to 1-electron term are periodic nearest-neighbor only");
    if (m_model_attrs.m_on_site_only_doubles)
        log::info("2-electron term diagonal contributions are on-site only");
}
