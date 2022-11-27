//
// Created by Robert J. Anderson on 09/01/2022.
//

#include "test_core/defs.h"
#include <test_core/sparse/Examples.h>
#include <M7_lib/linalg/Dense.h>
#include "M7_lib/arnoldi/ArnoldiSolver.h"
#include "M7_lib/util/Sort.h"
#include "M7_lib/hdf5/File.h"

TEST(ArnoldiSolver, SymNonDist) {
    const uint_t nrow = 20;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.symmetrized(false);

    ArnoldiOptions opts;
    opts.m_nroot = 3ul;

    ArnoldiSolver<double> solver(sym, nrow, opts, ArnoldiSolverBase::c_sym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    dense::SquareMatrix<double> dense(sym);
    v_t<double> evals;
    dense::diag(dense, evals);
    // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
    sort::inplace(evals, false, true);
    double eval;
    for (uint_t iroot=0ul; iroot<opts.m_nroot; ++iroot){
        solver.get_eval(iroot, eval);
        ASSERT_NEAR_EQ(eval, evals[iroot]);
    }
}

TEST(ArnoldiSolver, SymDist) {
    const uint_t nrow = 20;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.symmetrized(false);
    dist_mv_prod::Sparse<double> prod(sym);

    ArnoldiOptions opts;
    opts.m_nroot = 3ul;

    ArnoldiSolver<double> solver(prod, opts, ArnoldiSolverBase::c_sym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    if (mpi::i_am_root()) {
        dense::SquareMatrix<double> dense(sym);
        v_t<double> evals;
        dense::diag(dense, evals);
        // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
        sort::inplace(evals, false, true);
        auto dense_eval_it = evals.cbegin();
        double eval;
        for (uint_t iroot = 0ul; iroot < opts.m_nroot; ++iroot) {
            solver.get_eval(iroot, eval);
            ASSERT_NEAR_EQ(eval, dense_eval_it[iroot]);
        }
    }
}

TEST(ArnoldiSolver, NonSymNonDist) {
    const uint_t nrow = 20;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);

    ArnoldiOptions opts;
    opts.m_nroot = 2ul;
    ArnoldiSolver<double> solver(mat, nrow, opts, ArnoldiSolverBase::c_nonsym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    dense::SquareMatrix<double> dense(mat);
    // complex eigenvalues required for non-symmetric matrix diagonalization
    v_t<std::complex<double>> evals;
    ASSERT_TRUE(dense::diag(dense, evals));
    // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
    sort::inplace(evals, false, true);
    auto dense_eval_it = evals.cbegin();
    std::complex<double> eval;
    for (uint_t iroot = 0ul; iroot < opts.m_nroot; ++iroot) {
        solver.get_eval(iroot, eval);
        ASSERT_NEAR_EQ(eval, dense_eval_it[iroot]);
    }
}

TEST(ArnoldiSolver, NonSymDist) {
    const uint_t nrow = 20;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    dist_mv_prod::Sparse<double> prod(mat);

    ArnoldiOptions opts;
    opts.m_nroot = 3ul;
    ArnoldiSolver<double> solver(prod, opts, ArnoldiSolverBase::c_nonsym);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    if (mpi::i_am_root()) {
        dense::SquareMatrix<double> dense(mat);
        v_t<std::complex<double>> evals;
        dense::diag(dense, evals);
        // sort the eigenvalues by magnitude, largest first, since this is the order found by ARPACK
        sort::inplace(evals, false, true);
        auto dense_eval_it = evals.cbegin();
        std::complex<double> eval;

        for (uint_t iroot = 0ul; iroot < opts.m_nroot; ++iroot) {
            solver.get_eval(iroot, eval);
            ASSERT_NEAR_EQ(eval, arith::real(dense_eval_it[iroot]));
        }
    }
}
