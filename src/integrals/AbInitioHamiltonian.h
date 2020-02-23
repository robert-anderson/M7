//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_ABINITIOHAMILTONIAN_H
#define M7_ABINITIOHAMILTONIAN_H


#include <cstddef>
#include <string>
#include "Integrals_1e.h"
#include "Integrals_2e.h"
#include "src/fermion/Determinant.h"

class AbInitioHamiltonian {
    FcidumpFileIterator<defs::ham_t> m_file_iterator;
    defs::ham_t m_int_0;
    Integrals_1e<defs::ham_t, defs::isym_1e> m_int_1;
    Integrals_2e<defs::ham_t, defs::isym_2e> m_int_2;

public:

    inline size_t norb() const { return m_file_iterator.m_norb; }

    inline size_t nspatorb() const { return m_int_1.m_nspatorb; }

    inline size_t nelec() const { return m_file_iterator.m_nelec; }

    inline size_t spin_resolved() const { return m_file_iterator.m_spin_resolved; }

    inline size_t spin_conserving() const { return m_int_1.spin_conserving(); }

    auto &int_1() const {return m_int_1;}
    auto &int_2() const {return m_int_2;}

    size_t m_nci;

    AbInitioHamiltonian(const std::string &fname);

    consts::component_t<defs::ham_t>::type get_energy(const Determinant &det) const;

    defs::ham_t get_element_0(const Determinant &det) const;

    defs::ham_t get_element_1(const Determinant &bra, const size_t &removed, const size_t &inserted) const;

    defs::ham_t get_element_1(const Determinant &bra, const Determinant &ket) const;

    defs::ham_t get_element_2(const Determinant &bra,
                              const size_t &removed1, const size_t &removed2,
                              const size_t &inserted1, const size_t &inserted2) const;

    defs::ham_t get_element_2(const Determinant &bra, const Determinant &ket) const;

    defs::ham_t get_element(const Determinant &bra, const Determinant &ket) const;

    Determinant guess_reference(const size_t &spin_level) const;

    Determinant refine_guess_reference(const Determinant ref) const;

    Determinant choose_reference(const size_t &spin_level) const;
};


#endif //M7_ABINITIOHAMILTONIAN_H
