//
// Created by rja on 27/02/2020.
//

#include <src/enumerators/BitfieldEnumerator.h>
#include <src/enumerators/VectorCombinationEnumerator.h>
#include "ExactPropagator.h"


ExactPropagator::ExactPropagator(const std::unique_ptr<Hamiltonian> &ham,
                                 const RankAllocator<Determinant> &rankAllocator, double tau, defs::ham_comp_t shift)
        : DeterministicPropagator(ham, rankAllocator, tau, shift) {}


void ExactPropagator::off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                                   const NumericView<bool> flag_deterministic, const NumericView<bool> flag_initiator,
                                   TableArray<SpawnList> &spawn_list) const {
    auto occs = DeterminantSetEnumerator(determinant).enumerate();
    assert(occs.size());
    auto unoccs = DeterminantClrEnumerator(determinant).enumerate();
    assert(unoccs.size());

    Determinant excited(determinant.nspatorb());
    for (auto occ : occs) {
        for (auto unocc :unoccs) {
            excited = determinant.get_excited_det(occ, unocc);
            auto delta = -*weight*m_tau*m_ham->get_element(determinant, excited);
            add_to_spawn_list(excited, delta, flag_initiator, spawn_list);
        }
    }

    VectorCombinationEnumerator occ_enumerator(occs, 2);
    defs::inds occ_inds(2);

    while (occ_enumerator.next(occ_inds)) {
        {
            VectorCombinationEnumerator unocc_enumerator(unoccs, 2);
            defs::inds unocc_inds(2);
            while (unocc_enumerator.next(unocc_inds)) {
                excited = determinant.get_excited_det(occ_inds, unocc_inds);
                auto delta = -*weight*m_tau*m_ham->get_element(determinant, excited);
                add_to_spawn_list(excited, delta, flag_initiator, spawn_list);
            }
        }
    }
}