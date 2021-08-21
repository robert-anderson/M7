//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "src/core/basis/DecodedDeterminant.h"
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

FermionHamiltonian::FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                                       bool spin_conserving_2e, bool complex_valued, bool spin_resolved, size_t int_2e_rank) :
        m_nelec(nelec), m_nsite(nsite),
        m_spin_conserving_1e(spin_conserving_1e),
        m_spin_conserving_2e(spin_conserving_2e),
        m_complex_valued(complex_valued),
        m_int_2e_rank(int_2e_rank),
        m_int_1(nsite, spin_resolved),
        m_int_2(nsite, spin_resolved)
{}

FermionHamiltonian::FermionHamiltonian(const FcidumpFileReader &file_reader) :
        FermionHamiltonian(file_reader.m_nelec, file_reader.m_nspatorb,
                           file_reader.m_spin_conserving_1e,
                           file_reader.m_spin_conserving_2e,
                           file_reader.m_complex_valued,
                           file_reader.m_spin_resolved,
                           file_reader.m_int_2e_rank) {

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;

    // assume no contribution cases to be nonzero unless a counterexample is found
    m_nonzero_contribs.fill(false);
    REQUIRE_LE_ALL(m_nelec, 2*m_nsite, "unphysical number of electrons specified in FCIDUMP file");

    log::info("Reading fermion Hamiltonian coefficients from FCIDUMP file \"" + file_reader.m_fname + "\"...");
    while (file_reader.next(inds, value)) {
        auto contrib_case = ham_data::get_coupling_contrib_case(inds);
        auto rank = c_contrib_inds[contrib_case].m_term.m_nann;
        // all integrals are particle number conserving (for now)
        ASSERT(rank == c_contrib_inds[contrib_case].m_term.m_ncre);
        m_nonzero_contribs[contrib_case] = true;
        if (rank==2) {
            if (m_on_site_only_0022 && !on_site(inds[0], inds[1], inds[2], inds[3]))
                m_on_site_only_0022 = false;
            m_int_2.set(inds, value);
        }
        else if (rank==1) {
            if (m_nn_only_1111 && !nearest_neighbors(inds[0], inds[1], false)) m_nn_only_1111 = false;
            if (m_nnp_only_1111 && !nearest_neighbors(inds[0], inds[1], true)) m_nnp_only_1111 = false;
            m_int_1.set(inds, value);
        }
        else if (rank==0) m_int_0 = value;
        else MPI_ABORT("File reader error");
    }
    mpi::barrier();
    log::info("FCIDUMP loading complete.");

    // some useful output to identify the kind of H detected
    if (!m_nonzero_contribs[contrib_0011])
        log::info("1-electron term in Hamiltonian has no diagonal contributions");
    if (!m_nonzero_contribs[contrib_1111])
        log::info("1-electron term in Hamiltonian has no single-excitation contributions");
    if (!m_nonzero_contribs[contrib_0022])
        log::info("2-electron term in Hamiltonian has no diagonal contributions");
    if (!m_nonzero_contribs[contrib_1122])
        log::info("2-electron term in Hamiltonian has no single-excitation contributions");
    if (!m_nonzero_contribs[contrib_2222])
        log::info("2-electron term in Hamiltonian has no double-excitation contributions");
    if (m_nnp_only_1111)
        log::info("single-excitation contributions to 1-electron term are periodic nearest-neighbor only");
    else if (m_nn_only_1111)
        log::info("single-excitation contributions to 1-electron term are nearest-neighbor only");
    if (m_on_site_only_0022)
        log::info("2-electron term diagonal contributions are on-site only");
}

FermionHamiltonian::FermionHamiltonian(std::string fname, bool spin_major) :
        FermionHamiltonian(FcidumpFileReader(fname, spin_major)){}


defs::ham_t FermionHamiltonian::get_element_0000(const field::FrmOnv &onv) const {
    defs::ham_t element = m_int_0;
    auto singles_fn = [&](const size_t& i){ element+=m_int_1(i, i);};
    auto doubles_fn = [&](const size_t& i, const size_t& j){ element+=m_int_2.phys_antisym_element(i, j, i, j);};
    onv.foreach_pair(singles_fn, doubles_fn);
    return element;
}


defs::ham_t FermionHamiltonian::get_element_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}