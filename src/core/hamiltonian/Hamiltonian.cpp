//
// Created by rja on 27/02/2020.
//

#include "Hamiltonian.h"

#if 0
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

Determinant Hamiltonian::choose_reference(const int &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}

MappedList<Determinant> Hamiltonian::all_connections_of_det(const Determinant &ref, const defs::ham_comp_t eps) const {
#if 0
    Specification spec;
    spec.add<Determinant>(m_nsite);
    spec.add<defs::ham_t>(1);
    class ConnectionList : public MappedList<Determinant> {
    public:
        Field <Determinant> m_determinant;
        Field <defs::ham_t> m_helement;

        ConnectionList(size_t nsite, size_t nbucket) : MappedList(nbucket, 0),
                                                       m_determinant(this, nsite), m_helement(this) {}
    };
    auto nbucket = integer_utils::combinatorial(nsite(), nelec());
    ConnectionList list(m_nsite, nbucket);

    auto occs = DeterminantSetEnumerator(ref).enumerate();
    auto unoccs = DeterminantClrEnumerator(ref).enumerate();

    Determinant excited(m_nsite);
    for (auto occ:occs) {
        for (auto unocc:unoccs) {
            excited = ref.get_excited_det(occ, unocc);
            auto helement = get_element(ref, excited);
            if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                size_t irow = list.push(excited);
                *list.view<defs::ham_t>(irow) = helement;
            }
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
                auto helement = get_element(ref, excited);
                if (!consts::float_nearly_zero(std::abs(helement), eps)) {
                    size_t irow = list.push(excited);
                    *list.view<defs::ham_t>(irow) = helement;
                    assert(list.lookup(list.view<Determinant>(irow, 0)) == irow);
                }
            }
        }
    }
    return list;
#endif
}
#endif
