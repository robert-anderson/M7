//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include "../defs.h"
#include "Matrix.h"
#include "../integrals/AbInitioHamiltonian.h"

class DenseHamiltonian : public Matrix<defs::ham_t, true>{
public:
    DenseHamiltonian(const AbInitioHamiltonian &source);
};


#endif //M7_DENSEHAMILTONIAN_H
