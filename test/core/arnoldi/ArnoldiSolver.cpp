//
// Created by rja on 09/01/2022.
//

#include <test/core/sparse/Examples.h>
#include "gtest/gtest.h"
#include "src/core/arnoldi/ArnoldiSolver.h"
#include "src/core/linalg/EigenSolver.h"

TEST(ArnoldiSolver, Test) {
    const size_t nrow = 7;
    auto mat = sparse_matrix_examples::rect_double(nrow, nrow, 2);
    auto sym = mat.get_symmetrized(false);
    dist_mv_prod::Sparse<double> prod(sym);

    ArnoldiProblemSym<double> arnoldi_problem(2);
    arnoldi_problem.solve(prod);

    EigenSolver<double> dense_solver(sym.to_dense());
    std::cout << dense_solver.m_evals << std::endl;
    std::cout << arnoldi_problem.real_eigenvalue(0) << std::endl;
    std::cout << arnoldi_problem.real_eigenvalue(1) << std::endl;
}