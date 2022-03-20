//
// Created by rja on 09/01/2022.
//

#include <test_core/sparse/Examples.h>
#include <M7_lib/linalg/Dense.h>
#include "M7_lib/arnoldi/ArnoldiSolver.h"
#include "gtest/gtest.h"

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
    dense::SquareMatrix<double> dense(sym);
    std::vector<double> evals;
    dense::diag(dense, evals);
    // Arnoldi finds the extremal eigenvalue. in this case, the most positive
    auto dense_eval_it = evals.cbegin()+(nrow-nroot);
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
        dense::SquareMatrix<double> dense(sym);
        std::vector<double> evals;
        dense::diag(dense, evals);
        auto dense_eval_it = evals.cbegin() + (nrow - nroot);
        for (size_t iroot = 0ul; iroot < nroot; ++iroot) {
            ASSERT_FLOAT_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
        }
    }
}
TEST(ArnoldiSolver, NonSymNonDist) {
    const size_t nrow = 20;
    const size_t nroot = 3;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);

    ArnoldiProblemNonSym<double> arnoldi_problem(nroot);
    arnoldi_problem.solve(mat);

    /*
     * check Arnoldi solution against dense LAPACK full diagonalization
     */
    //EigenSolver<double> dense_solver(mat.to_dense());
    std::cout << dense::SquareMatrix<double>(mat).to_string() << std::endl;
    //auto dense_eval_it = dense_solver.m_evals.cbegin()+(nrow-nroot);
    for (size_t iroot=0ul; iroot<nroot; ++iroot){
        std::cout << arnoldi_problem.real_eigenvalue(iroot) << std::endl;
        //ASSERT_FLOAT_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
    }
}

TEST(ArnoldiSolver, NonSymDist) {
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
        dense::SquareMatrix<double> dense(sym);
        std::vector<double> evals;
        dense::diag(dense, evals);
        auto dense_eval_it = evals.cbegin() + (nrow - nroot);
        for (size_t iroot = 0ul; iroot < nroot; ++iroot) {
            ASSERT_FLOAT_EQ(arnoldi_problem.real_eigenvalue(iroot), dense_eval_it[iroot]);
        }
    }
}
