//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include "src/core/basis/DeterminantList.h"
#include "src/core/util/defs.h"
#include "Matrix.h"
#include "src/core/hamiltonian/Hamiltonian.h"

class DenseHamiltonian : public Matrix<defs::ham_t> {
public:
    DenseHamiltonian(const Hamiltonian &source);
    DenseHamiltonian(const Hamiltonian &source, DeterminantList &detlist);
};

#endif //M7_DENSEHAMILTONIAN_H
