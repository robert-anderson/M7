//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include "src/defs.h"
#include "Matrix.h"
#include "src/hamiltonian/Hamiltonian.h"

class DenseHamiltonian : public Matrix<defs::ham_t>{
public:
    DenseHamiltonian(const Hamiltonian &source);
};


#endif //M7_DENSEHAMILTONIAN_H
