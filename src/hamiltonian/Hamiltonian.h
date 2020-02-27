//
// Created by rja on 27/02/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <cstddef>
#include "src/consts.h"
#include "src/defs.h"
#include "src/fermion/Determinant.h"

class Hamiltonian {
    const size_t m_nsite;

public:
    Hamiltonian(const size_t &nsite);

    virtual consts::component_t<defs::ham_t>::type get_energy(const Determinant &det) const = 0;

    virtual defs::ham_t get_element_0(const Determinant &det) const = 0;

    virtual defs::ham_t get_element_1(const Determinant &bra, const size_t &removed, const size_t &inserted) const = 0;

    virtual defs::ham_t get_element_1(const Determinant &bra, const Determinant &ket) const = 0;

    virtual defs::ham_t get_element_2(const size_t &removed1, const size_t &removed2,
                              const size_t &inserted1, const size_t &inserted2) const = 0;

    virtual defs::ham_t get_element_2(const Determinant &bra, const Determinant &ket) const = 0;

    virtual defs::ham_t get_element(const Determinant &bra, const Determinant &ket) const = 0;

    size_t nsite() const {
        return m_nsite;
    }

    virtual bool spin_conserving() const = 0;

};


#endif //M7_HAMILTONIAN_H
