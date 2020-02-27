//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/enumerators/VectorCombinationEnumerator.h>
#include "AbInitioHamiltonian.h"
#include "../io/FcidumpFileIterator.h"
#include "../enumerators/BitfieldEnumerator.h"

AbInitioHamiltonian::AbInitioHamiltonian(const std::string &fname) :
    Hamiltonian(m_file_iterator.m_norb),
    m_file_iterator(fname),
    m_int_1{m_file_iterator.m_norb, spin_resolved()},
    m_int_2{m_file_iterator.m_norb, spin_resolved()} {
    defs::inds inds(4);
    defs::ham_t value;
    while (m_file_iterator.next(inds, value)) {
        if (m_int_2.valid_inds(inds)) m_int_2.set_from_fcidump(inds, value);
        else if (m_int_1.valid_inds(inds)) m_int_1.set_from_fcidump(inds, value);
        else if (inds[0] == ((size_t) -1)) m_int_0 = value;
    }
    m_nci = integer_utils::combinatorial(norb(), nelec());
}

defs::ham_t AbInitioHamiltonian::get_element_0(const Determinant &det) const {
    defs::ham_t element = m_int_0;
    DeterminantSetEnumerator set_inds_outer(det);
    size_t i;
    while (set_inds_outer.next(i)) {
        element += m_int_1.get(i, i);
        {
            size_t j;
            DeterminantSetEnumerator set_inds_inner(det);
            while (set_inds_inner.next(j) && j < i) {
                element += m_int_2.get_phys_antisym(i, j, i, j);
            }
        }
    }
    return element;
}

defs::ham_t AbInitioHamiltonian::get_element_1(const Determinant &ket,
                                               const size_t &removed, const size_t &inserted) const {
    DeterminantSetEnumerator common_inds(ket);
    size_t common;
    defs::ham_t element = m_int_1.get(inserted, removed);
    while (common_inds.next(common)) {
        if (common == removed) continue;
        element += m_int_2.get_phys_antisym(inserted, common, removed, common);
    }
    return element;
}

defs::ham_t AbInitioHamiltonian::get_element_1(const Determinant &bra, const Determinant &ket) const {
    size_t removed, inserted;
    {
        DeterminantAndNotEnumerator enumerator(ket, bra);
        enumerator.next(removed);
    }
    {
        DeterminantAndNotEnumerator enumerator(bra, ket);
        enumerator.next(inserted);
    }
    return get_element_1(ket, removed, inserted);
}

defs::ham_t AbInitioHamiltonian::get_element_2(const Determinant &bra, const Determinant &ket) const {
    size_t removed1, removed2, inserted1, inserted2;
    {
        DeterminantAndNotEnumerator enumerator(ket, bra);
        enumerator.next(removed1);
        enumerator.next(removed2);
    }
    {
        DeterminantAndNotEnumerator enumerator(bra, ket);
        enumerator.next(inserted1);
        enumerator.next(inserted2);
    }
    return m_int_2.get_phys_antisym(inserted1, inserted2, removed1, removed2);
}

defs::ham_t AbInitioHamiltonian::get_element_2(
    const size_t &removed1, const size_t &removed2,
    const size_t &inserted1, const size_t &inserted2) const {
    return m_int_2.get_phys_antisym(inserted1, inserted2, removed1, removed2);
}


defs::ham_t AbInitioHamiltonian::get_element(const Determinant &bra, const Determinant &ket) const {
    bool phase;
    switch (bra.nexcit(ket)) {
        case 0:
            return get_element_0(bra);
        case 1:
            phase = bra.phase(ket);
            assert(phase == ket.phase(bra));
            return phase ? -get_element_1(bra, ket) : get_element_1(bra, ket);
        case 2:
            phase = bra.phase(ket);
            assert(phase == ket.phase(bra));
            return phase ? -get_element_2(bra, ket) : get_element_2(bra, ket);
        default:
            return 0;
    }
}

defs::ham_comp_t AbInitioHamiltonian::get_energy(const Determinant &det) const {
    return consts::real(get_element_0(det));
}

Determinant AbInitioHamiltonian::guess_reference(const size_t &spin_level) const {
    Determinant ref(nspatorb());
    for (size_t i = 0ul; i < nelec() / 2 + 2 * spin_level + nelec() % 2; ++i) ref.set(i, 0);
    for (size_t i = 0ul; i < nelec() / 2; ++i) ref.set(i, 1);
    return ref;
}

Determinant AbInitioHamiltonian::refine_guess_reference(const Determinant ref) const {

    auto e_ref = get_energy(ref);
    /*
    * check that none of the single and double connections have a lower energy
    */

    auto occs = DeterminantSetEnumerator(ref).enumerate();
    auto unoccs = DeterminantClrEnumerator(ref).enumerate();

    Determinant excited(nspatorb());
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

Determinant AbInitioHamiltonian::choose_reference(const size_t &spin_level) const {
    auto ref = guess_reference(spin_level);
    ref = refine_guess_reference(ref);
    return ref;
}