//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/core/util/defs.h>
#include "BosonCouplings.h"

using Hamiltonian = std::conditional<defs::bosons, BosonCouplings, FermionHamiltonian>::type;

#endif //M7_HAMILTONIAN_H
