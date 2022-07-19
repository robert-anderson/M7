//
// Created by anderson on 18/07/2022.
//

#ifndef M7_FCIINITIALIZER_H
#define M7_FCIINITIALIZER_H

#include <M7_lib/linalg/FciIters.h>
#include <M7_lib/arnoldi/ArnoldiSolver.h>
#include "Wavefunction.h"

struct FciInitializer {
    double m_eval;

    explicit FciInitializer(const Hamiltonian& h);
};


#endif //M7_FCIINITIALIZER_H
