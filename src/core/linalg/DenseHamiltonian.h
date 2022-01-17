//
// Created by Robert John Anderson on 2020-01-18.
//

#ifndef M7_DENSEHAMILTONIAN_H
#define M7_DENSEHAMILTONIAN_H

#include "src/defs.h"
#include "Dense.h"
#include "src/core/hamiltonian/Hamiltonian.h"

class DenseHamiltonian : public dense::SquareMatrix<defs::ham_t> {

    void setup_frm(const Hamiltonian &source);

    void setup_frmbos(const Hamiltonian &source);

    void setup_bos(const Hamiltonian &source);

public:

    DenseHamiltonian(const Hamiltonian &source);

};

#endif //M7_DENSEHAMILTONIAN_H
