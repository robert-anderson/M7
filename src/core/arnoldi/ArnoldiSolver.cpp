//
// Created by rja on 09/01/2022.
//

#include "ArnoldiSolver.h"

bool ArnoldiProblemBase::solve_base(const std::function<void()> &product_fn) {
    size_t nmv_call = 0ul;
    while (!basis_found()){
        take_step();
        if (do_another_mv_call()){
            ++nmv_call;
            product_fn();
        }
    }
    // Finding eigenvalues and eigenvectors.
    find_eigenvalues();
    log::info("{} converged after {}",
              string_utils::plural("eigenpair", m_nroot),
              string_utils::plural("parallel matrix-vector multiplication", nmv_call));
    return true;
}
