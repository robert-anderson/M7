//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include "src/defs.h"
#include "Matrix.h"
#include "src/core/hamiltonian/FermiBosHamiltonian.h"

class DenseHamiltonian : public Matrix<defs::ham_t> {
public:
    DenseHamiltonian(const FermionHamiltonian &source);
    DenseHamiltonian(const FermiBosHamiltonian &source, int spin);
    //DenseHamiltonian(const FermionHamiltonian &source, const BosonCouplings& bc);
    //DenseHamiltonian(const FermionHamiltonian &source, DeterminantList &detlist);
};

#endif //M7_DENSEHAMILTONIAN_H
