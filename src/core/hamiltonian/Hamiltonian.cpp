//
// Created by rja on 27/02/2020.
//

#include <src/core/enumerator/ContainerCombinationEnumerator.h>
#include <src/core/fermion/DecodedDeterminant.h>
#include <src/core/parallel/RankAllocator.h>
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(const size_t &nsite) : m_nsite(nsite) {

}

Determinant Hamiltonian::guess_reference(const int &spin_restrict) const {
    Determinant ref(m_nsite);
    ASSERT((size_t)abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set(0, i);
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set(1, i);
    ASSERT(ref.spin() == spin_restrict);
    return ref;
}

Determinant Hamiltonian::refine_guess_reference(const DeterminantElement &ref) const {

    auto e_ref = get_energy(ref);
    /*
    * check that none of the single and double connections have a lower energy
    */

    OccupiedOrbitals occs(ref);
    VacantOrbitals vacs(ref);

    Determinant excited(m_nsite);

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

Determinant Hamiltonian::choose_reference(const int &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}

void
Hamiltonian::all_connections_of_det(ConnectionList *list, const Determinant &ref, const defs::ham_comp_t eps) const {
    OccupiedOrbitals occs(ref);
    VacantOrbitals vacs(ref);
    AntisymConnection connection(ref);

    Determinant excited(m_nsite);
    for (size_t iocc = 0ul; iocc < occs.m_nind; ++iocc) {
        const auto &occ = occs.m_inds[iocc];
        for (size_t ivac = 0ul; ivac < vacs.m_nind; ++ivac) {
            const auto &vac = vacs.m_inds[ivac];
            connection.zero();
            connection.add(occ, vac);
            connection.apply(ref, excited);
            auto helement = get_element(connection);
            if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                size_t irow = list->push(excited);
                list->helement(irow) = helement;
            }
        }
    }

    ContainerCombinationEnumerator<defs::det_work> occ_enumerator(occs.m_inds, occs.m_nind, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            ContainerCombinationEnumerator<defs::det_work> vac_enumerator(vacs.m_inds, vacs.m_nind, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {
                connection.zero();
                connection.add(occ_inds[0], occ_inds[1], vac_inds[0], vac_inds[1]);
                connection.apply(ref, excited);
                auto helement = get_element(connection);
                if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                    size_t irow = list->push(excited);
                    list->helement(irow) = helement;
                    ASSERT(list->lookup(list->determinant(irow)) == irow);
                }
            }
        }
    }
}

void
Hamiltonian::generate_ci_space(WalkerList *list, RankAllocator<DeterminantElement> &ra, const int &spin_level) const {
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

    Determinant work_det(nsite());
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
