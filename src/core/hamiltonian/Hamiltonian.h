//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/defs.h>
#include "FermiBosHamiltonian.h"

template <bool enable_bosons = defs::enable_bosons>
using Hamiltonian = typename std::conditional<enable_bosons, FermiBosHamiltonian, FermionHamiltonian>::type;

#endif //M7_HAMILTONIAN_H
