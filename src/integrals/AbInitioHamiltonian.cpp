//
// Created by Robert John Anderson on 2020-01-18.
//

#include "AbInitioHamiltonian.h"
#include "../io/FcidumpFileIterator.h"
#include "../enumerators/BitfieldEnumerator.h"

AbInitioHamiltonian::AbInitioHamiltonian(const std::string &fname) :
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
    while(common_inds.next(common)){
        if (common==removed) continue;
        element+=m_int_2.get_phys_antisym(inserted, common, removed, common);
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

defs::ham_t AbInitioHamiltonian::get_element(const Determinant &bra, const Determinant &ket) const {
    bool phase;
    switch (bra.nexcit(ket)){
        case 0: return get_element_0(bra);
        case 1:
            phase = bra.phase(ket);
            assert(phase==ket.phase(bra));
            return phase?-get_element_1(bra, ket):get_element_1(bra, ket);
        case 2:
            phase = bra.phase(ket);
            assert(phase==ket.phase(bra));
            return phase?-get_element_2(bra, ket):get_element_2(bra, ket);
        default: return 0;
    }
}
