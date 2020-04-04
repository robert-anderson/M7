//
// Created by Robert John Anderson on 2020-01-18.
//

#include <src/core/fermion/DecodedDeterminant.h>
#include <src/core/fermion/Connection.h>
#include "AbInitioHamiltonian.h"
#include "src/core/io/FcidumpFileIterator.h"

AbInitioHamiltonian::AbInitioHamiltonian(const std::string &fname) :
    Hamiltonian(FcidumpFileIterator<defs::ham_t>(fname).nsite()),
    m_file_iterator(fname),
    m_int_1(m_file_iterator.m_norb, spin_resolved()),
    m_int_2(m_file_iterator.m_norb, spin_resolved()) {
    defs::inds inds(4);
    defs::ham_t value;
    while (m_file_iterator.next(inds, value)) {
        if (m_int_2.valid_inds(inds)) m_int_2.set_from_fcidump(inds, value);
        else if (m_int_1.valid_inds(inds)) m_int_1.set_from_fcidump(inds, value);
        else if (inds[0] == ((size_t) -1)) m_int_0 = value;
    }
    m_nci = integer_utils::combinatorial(nsite() * 2, nelec());
}

defs::ham_t AbInitioHamiltonian::get_element_0(const DeterminantElement &det) const {
    defs::ham_t element = m_int_0;
    OccupiedOrbitals occs(det);
    for (size_t i=0ul; i<occs.m_nind; ++i) {
        const auto &icommon = occs.m_inds[i];
        element += m_int_1.element(icommon, icommon);
        for (size_t j=0ul; j<i; ++j) {
            const auto &jcommon = occs.m_inds[j];
            element += m_int_2.phys_antisym_element(icommon, jcommon, icommon, jcommon);
        }
    }
    return element;
}

defs::ham_t AbInitioHamiltonian::get_element_1(const DeterminantElement &ket,
                                               const size_t &removed, const size_t &inserted) const {
    OccupiedOrbitals occs(ket);
    defs::ham_t element = m_int_1.element(inserted, removed);

    for (size_t i=0ul; i<occs.m_nind; ++i) {
        const auto &icommon = occs.m_inds[i];
        if (icommon == removed) continue;
        assert(icommon!=removed);
        assert(icommon!=inserted);
        element += m_int_2.phys_antisym_element(inserted, icommon, removed, icommon);
    }
    return element;
}

defs::ham_t AbInitioHamiltonian::get_element_2(
    const size_t &removed1, const size_t &removed2,
    const size_t &inserted1, const size_t &inserted2) const {
    return m_int_2.phys_antisym_element(inserted1, inserted2, removed1, removed2);
}

defs::ham_t AbInitioHamiltonian::get_element(
    const DeterminantElement &ket, const AntisymConnection &connection) const {
    defs::ham_t helement = 0;
    const auto &cre = connection.m_cre;
    const auto &des = connection.m_des;
    switch (connection.nexcit()) {
        case 0:
            return get_element_0(ket);
        case 1:
            helement = get_element_1(ket, des[0], cre[0]);
            break;
        case 2:
            helement = get_element_2(des[0], des[1], cre[0], cre[1]);
        default:
            return 0;
    }
    return connection.m_phase ? -helement : helement;
}

defs::ham_comp_t AbInitioHamiltonian::get_energy(const DeterminantElement &det) const {
    return consts::real(get_element_0(det));
}

size_t AbInitioHamiltonian::nelec() const { return m_file_iterator.m_nelec; }

bool AbInitioHamiltonian::spin_resolved() const { return m_file_iterator.m_spin_resolved; }

bool AbInitioHamiltonian::spin_conserving() const { return m_int_1.spin_conserving(); }

auto &AbInitioHamiltonian::int_1() const { return m_int_1; }

auto &AbInitioHamiltonian::int_2() const { return m_int_2; }