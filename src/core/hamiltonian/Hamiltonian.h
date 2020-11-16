//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/core/util/defs.h>
#include "FermiBosHamiltonian.h"

using Hamiltonian = std::conditional<defs::bosons, FermiBosHamiltonian, FermionHamiltonian>::type;

#endif //M7_HAMILTONIAN_H
