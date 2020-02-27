//
// Created by rja on 27/02/2020.
//

#include <src/enumerators/BitfieldEnumerator.h>
#include <src/enumerators/VectorCombinationEnumerator.h>
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(const size_t &nsite) : m_nsite(nsite){

}

Determinant Hamiltonian::guess_reference(const size_t &spin_level) const {
    Determinant ref(m_nsite);
    for (size_t i = 0ul; i < nelec() / 2 + 2 * spin_level + nelec() % 2; ++i) ref.set(i, 0);
    for (size_t i = 0ul; i < nelec() / 2; ++i) ref.set(i, 1);
    return ref;
}

Determinant Hamiltonian::refine_guess_reference(const Determinant ref) const {

    auto e_ref = get_energy(ref);
    /*
    * check that none of the single and double connections have a lower energy
    */

    auto occs = DeterminantSetEnumerator(ref).enumerate();
    auto unoccs = DeterminantClrEnumerator(ref).enumerate();

    Determinant excited(m_nsite);
    for (auto occ:occs) {
        for (auto unocc:unoccs) {
            excited = ref.get_excited_det(occ, unocc);
            if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
        }
    }

    VectorCombinationEnumerator occ_enumerator(occs, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            VectorCombinationEnumerator unocc_enumerator(unoccs, 2);
            defs::inds unocc_inds(2);
            while (unocc_enumerator.next(unocc_inds)) {
                excited = ref.get_excited_det(occ_inds, unocc_inds);
                if (get_energy(excited) < e_ref) return refine_guess_reference(excited);
            }
        }
    }
    return ref;
}

Determinant Hamiltonian::choose_reference(const size_t &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}