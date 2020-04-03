//
// Created by rja on 27/02/2020.
//

#include <src/core/enumerator/VectorCombinationEnumerator.h>
#include <src/core/fermion/DecodedDeterminant.h>
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(const size_t &nsite) : m_nsite(nsite) {

}

Determinant Hamiltonian::guess_reference(const int &spin_restrict) const {
    Determinant ref(m_nsite);
    assert(abs(spin_restrict) % 2 == nelec() % 2);
    size_t n_spin_0 = (nelec() + spin_restrict) / 2;
    size_t n_spin_1 = nelec() - n_spin_0;
    for (size_t i = 0ul; i < n_spin_0; ++i) ref.set(i, 0);
    for (size_t i = 0ul; i < n_spin_1; ++i) ref.set(i, 1);
    assert(ref.spin() == spin_restrict);
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

    VectorCombinationEnumerator occ_enumerator(occs.m_inds, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            VectorCombinationEnumerator vac_enumerator(vacs.m_inds, 2);
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

Hamiltonian::ConnectionList
Hamiltonian::all_connections_of_det(const Determinant &ref, const defs::ham_comp_t eps) const {

    auto nbucket = integer_utils::combinatorial(nsite(), nelec());
    ConnectionList list(m_nsite, nbucket);

    OccupiedOrbitals occs(ref);
    VacantOrbitals vacs(ref);

    Determinant excited(m_nsite);
    for (auto occ : occs.m_inds) {
        for (auto vac : vacs.m_inds) {
            excited = ref;
            excited.excite(vac, occ);
            auto helement = get_element(ref, excited);
            if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                size_t irow = list.push(excited);
                list.helement.element(irow) = helement;
            }
        }
    }

    VectorCombinationEnumerator occ_enumerator(occs.m_inds, 2);
    defs::inds occ_inds(2);
    while (occ_enumerator.next(occ_inds)) {
        {
            VectorCombinationEnumerator vac_enumerator(vacs.m_inds, 2);
            defs::inds vac_inds(2);
            while (vac_enumerator.next(vac_inds)) {
                excited = ref;
                excited.excite(vac_inds[0], vac_inds[1], occ_inds[0], occ_inds[1]);
                auto helement = get_element(ref, excited);
                if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                    size_t irow = list.push(excited);
                    list.helement.element(irow) = helement;
                    assert(list.lookup(list.helement.element(irow)) == irow);
                }
            }
        }
    }
    return list;
}
