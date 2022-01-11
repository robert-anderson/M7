//
// Created by rja on 09/01/2022.
//

#include <test/core/sparse/Examples.h>
#include "gtest/gtest.h"
#include "src/core/arnoldi/ArnoldiSolver.h"
#include "src/core/linalg/EigenSolver.h"

TEST(ArnoldiSolver, SymNonDist) {
    const size_t nrow = 20;
    const size_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.get_symmetrized(false);

    ArnoldiProblemSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(sym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    EigenSolver<double> dense_solver(sym.to_dense());
    auto dense_eval_it = dense_solver.m_evals.cbegin()+(nrow-nroot);
    for (size_t iroot=0ul; iroot<nroot; ++iroot){
        ASSERT_FLOAT_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
    }
}

TEST(ArnoldiSolver, SymDist) {
    const size_t nrow = 20;
    const size_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.get_symmetrized(false);
    dist_mv_prod::Sparse<double> prod(sym);

    ArnoldiProblemSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(prod);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    if (mpi::i_am_root()) {
        EigenSolver<double> dense_solver(sym.to_dense());
        auto dense_eval_it = dense_solver.m_evals.cbegin() + (nrow - nroot);
        for (size_t iroot = 0ul; iroot < nroot; ++iroot) {
            ASSERT_FLOAT_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
        }
    }
}