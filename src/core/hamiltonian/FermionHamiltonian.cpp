//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "src/core/basis/DecodedDeterminant.h"
#include "src/core/parallel/RankAllocator.h"
#include "FermionHamiltonian.h"


buffered::FermionOnv FermionHamiltonian::guess_reference(const int &spin_restrict) const {
    buffered::FermionOnv ref(m_nsite);
    ASSERT((size_t)abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set({0, i});
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set({1, i});
    ASSERT(ref.spin() == spin_restrict);
    return ref;
}
#if 0

FermionOnv FermionHamiltonian::refine_guess_reference(const DeterminantElement &ref) const {

    auto e_ref = get_energy(ref);
    /*
    * check that none of the single and double connections have a lower energy
    */

    OccupiedOrbitals occs(ref);
    VacantOrbitals vacs(ref);

    FermionOnv excited(m_nsite);

    for (auto occ : occs.m_inds) {
        for (auto vac : vacs.m_inds) {
            excited = ref;
            excited.excite(vac, occ);
            if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
        }
    }

    ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {
                excited = ref;
                excited.excite(vac_inds[0], vac_inds[1], occ_inds[0], occ_inds[1]);
                if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
            }
        }
    }
    excited = ref;
    return excited;
}

FermionOnv FermionHamiltonian::choose_reference(const int &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}

void
FermionHamiltonian::generate_ci_space(WalkerTable *list, RankAllocator<DeterminantElement> &ra, const int &spin_level) const {
    size_t nalpha = nelec() / 2 + spin_level;
    size_t nbeta = nelec() / 2 - spin_level;

    defs::inds alpha_sites(nsite(), 0);
    defs::inds beta_sites(nsite(), 0);

    for (size_t i = 0; i < nsite(); ++i) {
        alpha_sites[i] = i;
        beta_sites[i] = nsite() + i;
    }
    ContainerCombinationEnumerator<defs::inds> alpha_enumerator(alpha_sites, nsite(), nalpha);
    defs::inds alpha_inds(nalpha);
    defs::inds beta_inds(nbeta);

    FermionOnv work_det(nsite());
    while (alpha_enumerator.next(alpha_inds)) {
        ContainerCombinationEnumerator<defs::inds> beta_enumerator(beta_sites, nsite(), nbeta);
        while (beta_enumerator.next(beta_inds)) {
            work_det.zero();
            work_det.set(alpha_inds);
            work_det.set(beta_inds);
            if (!mpi::i_am(ra.get_rank(work_det))) continue;
            auto irow = list->expand_push(work_det);
            auto det = list->m_determinant(irow);
            auto h_diag = list->m_hdiag(irow);
            det.set(alpha_inds);
            det.set(beta_inds);
            h_diag = get_energy(det);
        }
    }
#ifndef NDEBUG
    if (mpi::nrank()==1) {
        ASSERT(list->high_water_mark(0) ==
               integer_utils::combinatorial(nsite(), nalpha) * integer_utils::combinatorial(nsite(), nbeta))
    }
#endif
}

#endif

FermionHamiltonian::FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                                       bool spin_conserving_2e, bool complex_valued, bool spin_resolved, size_t int_2e_rank) :
        m_nelec(nelec), m_nsite(nsite),
        m_spin_conserving_1e(spin_conserving_1e),
        m_spin_conserving_2e(spin_conserving_2e),
        m_complex_valued(complex_valued),
        m_int_2e_rank(int_2e_rank),
        m_int_1(nsite, spin_resolved),
        m_int_2(nsite, spin_resolved),
        m_sym_helper(new ham_sym_helpers::Fermion(*this))
{}

FermionHamiltonian::FermionHamiltonian(const FcidumpFileReader &file_reader) :
        FermionHamiltonian(file_reader.nelec(), file_reader.nspatorb(),
                           file_reader.spin_conserving_1e(),
                           file_reader.spin_conserving_2e(),
                           file_reader.m_complex_valued,
                           file_reader.spin_resolved(),
                           file_reader.int_2e_rank()) {

    using namespace ham_data;
    defs::inds inds(4);
    defs::ham_t value;

    // assume no contribution cases to be nonzero unless a counterexample is found
    m_nonzero_contribs.fill(false);

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

defs::ham_t FermionHamiltonian::get_element_0(const defs::inds &occs, const size_t &nocc) const {
    defs::ham_t element = m_int_0;
    for (size_t i = 0ul; i < nocc; ++i) {
        auto const &occi = occs[i];
        element += m_int_1(occi, occi);
        for (size_t j = 0ul; j < i; ++j) {
            auto const &occj = occs[j];
            element += m_int_2.phys_antisym_element(occi, occj, occi, occj);
        }
    }
    return element;
}

defs::ham_t FermionHamiltonian::get_element_0(const OccupiedOrbitals &occs) const {
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t FermionHamiltonian::get_element_0(const fields::Onv<0> &fonv) const {
    // TODO: make occs member data
    OccupiedOrbitals occs(fonv);
    return get_element_0(occs.inds(), occs.size());
}

defs::ham_t FermionHamiltonian::get_element_0(const conn::Antisym<0> &connection) const {
    return get_element_0(connection.com(), connection.ncom());
}

defs::ham_t FermionHamiltonian::get_element_1(const conn::Antisym<0> &connection) const {
    const auto &cre = connection.cre(0);
    const auto &ann = connection.ann(0);
    const auto &coms = connection.com();
    const auto &ncom = connection.ncom();

    defs::ham_t element = m_int_1(cre, ann);
    for (size_t icom = 0ul; icom < ncom; ++icom)
        element += m_int_2.phys_antisym_element(cre, coms[icom], ann, coms[icom]);
    return connection.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
    return m_int_2.phys_antisym_element(i, j, k, l);
}

defs::ham_t FermionHamiltonian::get_element_2(const conn::Basic<0> &connection) const {
    return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
}

defs::ham_t FermionHamiltonian::get_element_2(const conn::Antisym<0> &connection) const {
    const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    return connection.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const {
    return get_element(AntisymFermionOnvConnection(ket, bra));
}

#if 0
FermionHamiltonian::SpinTerms::SpinTerms(const FermionHamiltonian &ham) : Terms(ham),
                                                                          m_spin_occ_work(ham.nsite()), m_spin_vac_work(ham.nsite()){
}

void FermionHamiltonian::SpinTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                       const FermionHamiltonian::Terms::body_fn_t &body, bool get_h,
                                                       bool h_nonzero_only, bool include_diagonal) const {
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    m_spin_vac_work.update(src_onv);
    if (include_diagonal) perform_diagonal(src_onv, body, get_h, h_nonzero_only);
    // spin a->a, aa->aa
    foreach_connection_subset(src_onv, m_spin_occ_work[0], m_spin_vac_work[0], body, get_h, h_nonzero_only);
    // spin b->b, bb->bb
    foreach_connection_subset(src_onv, m_spin_occ_work[1], m_spin_vac_work[1], body, get_h, h_nonzero_only);
    // spin ab->ab
    foreach_connection_subset_doubles(src_onv, m_spin_occ_work[0], m_spin_occ_work[1],
                                      m_spin_vac_work[0], m_spin_vac_work[1], body, get_h, h_nonzero_only);
}

FermionHamiltonian::Hubbard1DTerms::Hubbard1DTerms(const FermionHamiltonian &ham) : SpinTerms(ham){}

void FermionHamiltonian::Hubbard1DTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                            const FermionHamiltonian::Terms::body_fn_t &body,
                                                            bool get_h, bool h_nonzero_only, bool include_diagonal) const {
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    if (include_diagonal) perform_diagonal(src_onv, body, get_h, h_nonzero_only);
    // spin a
    for (auto& occ: m_spin_occ_work[0]) {
        // to the left
        if (occ > 0 && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.m_nsite && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
    }
    // spin b
    for (auto& occ: m_spin_occ_work[1]) {
        // to the left
        if (occ > src_onv.m_nsite && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.nbit() && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
    }
}

void FermionHamiltonian::Hubbard1DPbcTerms::foreach_connection(const fields::FermionOnv &src_onv,
                                                               const FermionHamiltonian::Terms::body_fn_t &body,
                                                               bool get_h, bool h_nonzero_only, bool include_diagonal) const {
    if (include_diagonal) perform_diagonal(src_onv, body, get_h, h_nonzero_only);
    m_helement_work = 0.0;
    m_spin_occ_work.update(src_onv);
    // spin a
    for (auto& occ: m_spin_occ_work[0]) {
        // to the left
        if (occ > 0 && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
            // to the left (wrap-around)
        else if (occ == 0 && !src_onv.get(src_onv.m_nsite))
            perform_single(src_onv, occ, src_onv.m_nsite-1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.m_nsite && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
            // to the right (wrap-around)
        else if (occ+1 == src_onv.m_nsite && !src_onv.get(0))
            perform_single(src_onv, occ, 0, body, get_h, h_nonzero_only);
    }
    // spin b
    for (auto& occ: m_spin_occ_work[1]) {
        // to the left
        if (occ > src_onv.m_nsite && !src_onv.get(occ - 1))
            perform_single(src_onv, occ, occ - 1, body, get_h, h_nonzero_only);
            // to the left (wrap-around)
        else if (occ == src_onv.m_nsite && !src_onv.get(src_onv.nbit()-1))
            perform_single(src_onv, occ, src_onv.nbit()-1, body, get_h, h_nonzero_only);
        // to the right
        if (occ+1<src_onv.nbit() && !src_onv.get(occ+1))
            perform_single(src_onv, occ, occ+1, body, get_h, h_nonzero_only);
            // to the right (wrap-around)
        else if (occ+1 == src_onv.nbit() && !src_onv.get(src_onv.m_nsite-1))
            perform_single(src_onv, occ, src_onv.m_nsite-1, body, get_h, h_nonzero_only);
    }
}
#endif