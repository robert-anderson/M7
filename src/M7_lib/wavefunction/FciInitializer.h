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

    struct MbfOrderRow : Row {
        field::Mbf m_mbf;
        MbfOrderRow(sys::Basis basis): m_mbf(this, basis, "key"){}

        field::Mbf &key_field(){return m_mbf;}
    };

    BufferedTable<MbfOrderRow, true> m_mbf_order_table;
    explicit FciInitializer(const Hamiltonian& h, ham_comp_t shift=0.0, ArnoldiOptions opts={});
};


#endif //M7_FCIINITIALIZER_H
