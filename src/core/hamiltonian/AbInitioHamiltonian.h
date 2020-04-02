//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_ABINITIOHAMILTONIAN_H
#define M7_ABINITIOHAMILTONIAN_H


#include <cstddef>
#include <string>
#include "Hamiltonian.h"
#include "src/integrals/Integrals_1e.h"
#include "src/integrals/Integrals_2e.h"

#if 0
#include "src/fermion/Determinant.h"

class AbInitioHamiltonian : public Hamiltonian {
    FcidumpFileIterator<defs::ham_t> m_file_iterator;
    defs::ham_t m_int_0;
    Integrals_1e<defs::ham_t, defs::isym_1e> m_int_1;
    Integrals_2e<defs::ham_t, defs::isym_2e> m_int_2;

public:

    explicit AbInitioHamiltonian(const std::string &fname);

    defs::ham_comp_t get_energy(const Determinant &det) const override;

    defs::ham_t get_element_0(const Determinant &det) const override;

    defs::ham_t get_element_1(const Determinant &bra, const size_t &removed, const size_t &inserted) const override;

    defs::ham_t get_element_1(const Determinant &bra, const Determinant &ket) const override;

    defs::ham_t get_element_2(const size_t &removed1, const size_t &removed2, const size_t &inserted1,
                              const size_t &inserted2) const override;

    defs::ham_t get_element_2(const Determinant &bra, const Determinant &ket) const override;

    defs::ham_t get_element(const Determinant &bra, const Determinant &ket) const override;

    size_t nelec() const override;

    bool spin_resolved() const override;

    bool spin_conserving() const override;

    auto &int_1() const;

    auto &int_2() const;
};


#endif //M7_ABINITIOHAMILTONIAN_H
#endif //M7_ABINITIOHAMILTONIAN_H
