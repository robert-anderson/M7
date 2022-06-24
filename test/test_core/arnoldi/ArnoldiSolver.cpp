//
// Created by Robert J. Anderson on 09/01/2022.
//

#include "test_core/defs.h"
#include <test_core/sparse/Examples.h>
#include <M7_lib/linalg/Dense.h>
#include "M7_lib/arnoldi/ArnoldiSolver.h"
#include "M7_lib/util/Sort.h"

TEST(ArnoldiSolver, SymNonDist) {
    const uint_t nrow = 20;
    const uint_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.get_symmetrized(false);

    ArnoldiProblemSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(sym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    dense::SquareMatrix<double> dense(sym);
    std::vector<double> evals;
    dense::diag(dense, evals);
    // Arnoldi finds the extremal eigenvalue. in this case, the most positive
    auto dense_eval_it = evals.cbegin()+(nrow-nroot);
    for (uint_t iroot=0ul; iroot<nroot; ++iroot){
        ASSERT_NEARLY_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
    }
}

TEST(ArnoldiSolver, SymDist) {
    const uint_t nrow = 20;
    const uint_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.get_symmetrized(false);
    dist_mv_prod::Sparse<double> prod(sym);

    ArnoldiProblemSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(prod);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    if (mpi::i_am_root()) {
        dense::SquareMatrix<double> dense(sym);
        std::vector<double> evals;
        dense::diag(dense, evals);
        auto dense_eval_it = evals.cbegin() + (nrow - nroot);
        for (uint_t iroot = 0ul; iroot < nroot; ++iroot) {
            ASSERT_NEARLY_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
        }
    }
}

TEST(ArnoldiSolver, NonSymNonDist) {
    const uint_t nrow = 20;
    const uint_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);

    ArnoldiProblemNonSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(mat);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    dense::SquareMatrix<double> dense(mat);
    // complex eigenvalues required for non-symmetric matrix diagonalization
    std::vector<std::complex<double>> evals;
    ASSERT_TRUE(dense::diag(dense, evals));
    // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
    sort::inplace(evals, false, true);
    auto dense_eval_it = evals.cbegin() + (nrow - nroot);
    for (uint_t iroot = 0ul; iroot < nroot; ++iroot) {
        ASSERT_NEARLY_EQ(arnoldi_problem.complex_eigenvalue(iroot), dense_eval_it[iroot]);
    }
}

TEST(ArnoldiSolver, NonSymDist) {
    const uint_t nrow = 20;
    const uint_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    dist_mv_prod::Sparse<double> prod(mat);

    ArnoldiProblemNonSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(prod);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    if (mpi::i_am_root()) {
        dense::SquareMatrix<double> dense(mat);
        std::vector<std::complex<double>> evals;
        dense::diag(dense, evals);
        // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
        sort::inplace(evals, false, true);
        auto dense_eval_it = evals.cbegin() + (nrow - nroot);
        for (uint_t iroot = 0ul; iroot < nroot; ++iroot) {
            ASSERT_NEARLY_EQ(arnoldi_problem.real_eigenvalue(iroot), arith::real(dense_eval_it[iroot]));
        }
    }
}
