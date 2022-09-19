//
// Created by anderson on 18/07/2022.
//

#ifndef M7_FCIINITIALIZER_H
#define M7_FCIINITIALIZER_H

#include <M7_lib/linalg/FciIters.h>
#include <M7_lib/arnoldi/ArnoldiSolver.h>
#include "Wavefunction.h"


struct FciInitOptions : ArnoldiOptions {
    /**
     * shift to add to the diagonal elements of the sparse subspace Hamiltonian
     */
    ham_comp_t m_diag_shift = 0.0;
};

struct FciInitializer {
    struct MbfOrderRow : Row {
        field::Mbf m_mbf;
        MbfOrderRow(sys::Basis basis): m_mbf(this, basis, "key"){}

        field::Mbf &key_field(){return m_mbf;}
    };
    /**
     * mapped list of basis functions to aid in the setup of sparse H, and retain the physical meaning of its rows
     */
    BufferedTable<MbfOrderRow, true> m_mbf_order_table;
    /**
     * instance of the ARPACK wrapper which is only used if H is hermitian (Lanczos)
     */
    ArnoldiProblemSym<ham_t> m_arpack_sym;
    /**
     * instance of the ARPACK wrapper which is only used if H is non-hermitian
     */
    ArnoldiProblemNonSym<ham_t> m_arpack_nonsym;
    explicit FciInitializer(const Hamiltonian& h, FciInitOptions opts={});

    ArnoldiResults<ham_comp_t> get_results() {
        return m_arpack_sym.m_solver ? m_arpack_sym.get_results() : m_arpack_nonsym.get_results();
    }
};


#endif //M7_FCIINITIALIZER_H
