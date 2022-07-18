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

    FciInitializer(const Hamiltonian& h) {
        auto row_iters = FciIters::make(h);
        auto col_iters = FciIters::make(h);
        const auto count = row_iters.niter_single();

        auto count_local = mpi::evenly_shared_count(count);
        auto displ_local = mpi::evenly_shared_displ(count);

        buffered::Mbf row_mbf(h.m_basis);
        auto col_mbf = row_mbf;

        sparse::dynamic::Matrix<double> sparse_ham;
        sparse_ham.resize(count_local);


        uint_t irow = ~0ul;
        auto fn_row_loop = [&]() {
            ++irow;
            if ((irow < displ_local) || (irow >= (displ_local+count_local))) return;
            uint_t icol = ~0ul;
            auto fn_col_loop = [&]() {
                ++icol;
                auto helem = h.get_element(row_mbf, col_mbf);
                if (!ham::is_significant(helem)) return;
                sparse_ham.insert(irow-displ_local, {icol, helem});
            };
            col_iters.m_single->loop(col_mbf, fn_col_loop);
        };
        row_iters.m_single->loop(row_mbf, fn_row_loop);

//        hdf5::FileWriter fw("ham.h5");
//        dense::Matrix<double>(sparse_ham).save("data", fw);
        const uint_t nroot = 3;
        ArnoldiProblemNonSym<double> solver(nroot);
        solver.solve(sparse_ham);
        m_eval = solver.real_eigenvalue(nroot-1);
        logging::info("Non-symmetric arnoldi eigenvalue: {}", m_eval);
    }
};


#endif //M7_FCIINITIALIZER_H
