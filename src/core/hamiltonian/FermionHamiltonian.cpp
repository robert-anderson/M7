//
// Created by rja on 27/02/2020.
//

#include "src/core/enumerator/ContainerCombinationEnumerator.h"
#include "src/core/basis/DecodedDeterminant.h"
#include "src/core/parallel/RankAllocator.h"
#include "FermionHamiltonian.h"



#if 0
FermionOnv FermionHamiltonian::guess_reference(const int &spin_restrict) const {
    FermionOnv ref(m_nsite);
    ASSERT((size_t)abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set(0, i);
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set(1, i);
    ASSERT(ref.spin() == spin_restrict);
    return ref;
}

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
FermionHamiltonian::generate_ci_space(WalkerList *list, RankAllocator<DeterminantElement> &ra, const int &spin_level) const {
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
#ifndef DNDEBUG
    if (mpi::nrank()==1) {
        ASSERT(list->high_water_mark(0) ==
               integer_utils::combinatorial(nsite(), nalpha) * integer_utils::combinatorial(nsite(), nbeta))
    }
#endif
}

#endif

FermionHamiltonian::FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e,
                                       bool spin_conserving_2e, bool complex_valued, bool spin_resolved) :
        m_nelec(nelec), m_nsite(nsite),
        m_spin_conserving_1e(spin_conserving_1e),
        m_spin_conserving_2e(spin_conserving_2e),
        m_complex_valued(complex_valued),
        m_int_1(nsite, spin_resolved),
        m_int_2(nsite, spin_resolved)
{}

FermionHamiltonian::FermionHamiltonian(const FcidumpFileReader<defs::ham_t> &file_reader) :
        FermionHamiltonian(file_reader.nelec(), file_reader.nspatorb(),
                           file_reader.spin_conserving_1e(),
                           file_reader.spin_conserving_2e(),
                           file_reader.m_complex_valued, file_reader.spin_resolved()) {
    defs::inds inds(4);
    defs::ham_t value;

    logger::write("Loading ab-initio FermionHamiltonian from FCIDUMP...");
    while (file_reader.next(inds, value)) {
        if (ints2_t::valid_inds(inds)) m_int_2.set(inds, value);
        else if (ints1_t::valid_inds(inds)) m_int_1.set(inds, value);
        else if (inds[0] == ~0ul) m_int_0 = value;
    }
    mpi::barrier();
    logger::write("FCIDUMP loading complete.");
}

FermionHamiltonian::FermionHamiltonian(const std::string &fname, bool spin_major) :
        FermionHamiltonian(FcidumpFileReader<defs::ham_t>(fname, spin_major)){}

consts::component_t<defs::ham_t>::type FermionHamiltonian::get_energy(const views::FermionOnv &det) const {
    return consts::real(get_element_0(det));
}

defs::ham_t FermionHamiltonian::get_element_0(const defs::det_work &occs, const size_t &nocc) const {
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
    return get_element_0(occs.m_inds, occs.m_nind);
}

defs::ham_t FermionHamiltonian::get_element_0(const views::FermionOnv &det) const {
    OccupiedOrbitals occs(det);
    return get_element_0(occs.m_inds, occs.m_nind);
}

defs::ham_t FermionHamiltonian::get_element_0(const AntisymConnection &connection) const {
    return get_element_0(connection.com(), connection.ncom());
}

defs::ham_t FermionHamiltonian::get_element_1(const AntisymConnection &connection) const {
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

defs::ham_t FermionHamiltonian::get_element_2(const DeterminantConnection &connection) const {
    return get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
}

defs::ham_t FermionHamiltonian::get_element_2(const AntisymConnection &connection) const {
    const auto element = get_element_2(connection.cre(0), connection.cre(1), connection.ann(0), connection.ann(1));
    return connection.phase() ? -element : element;
}

defs::ham_t FermionHamiltonian::get_element(const AntisymConnection &connection) const {
    switch (connection.nexcit()) {
        case 0:
            return get_element_0(connection);
        case 1: ASSERT(connection.ncom() + connection.nexcit() == nelec());
            return get_element_1(connection);
        case 2:
            return get_element_2(connection);
        default:
            return 0;
    }
}

defs::ham_t FermionHamiltonian::get_element(const views::FermionOnv &bra, const views::FermionOnv &ket) const {
    return get_element(AntisymConnection(ket, bra));
}
