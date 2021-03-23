//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/defs.h>
#include "FermiBosHamiltonian.h"
#include "src/core/nd/NdArray.h"

template <bool enable_bosons = defs::enable_bosons>
using Hamiltonian = typename std::conditional<enable_bosons, FermiBosHamiltonian, FermionHamiltonian>::type;

/*
struct ConnectionGeneral {
    size_t m_cre_f;
    size_t m_ann_f;


};

struct HamiltonianGeneral {
    const std::vector<std::array<size_t, 4>> m_max_orders;
    // 4 operator types + 2 excit levels
    NdArray<bool, 6> m_nonzero_orders;
    const std::vector<defs::inds> m_exlvls;
    //HamiltonianGeneral(std::vector<std::array<size_t, 4>>& m_orders)
    HamiltonianGeneral(){
        m_nonzero_orders.m_format
    }

};
*/

#endif //M7_HAMILTONIAN_H
