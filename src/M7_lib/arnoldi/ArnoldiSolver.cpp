//
// Created by Robert J. Anderson on 09/01/2022.
//

#include "ArnoldiSolver.h"

bool ArnoldiProblemBase::solve_base(const std::function<void()> &product_fn, bool dist) {
    bool i_am_solver_rank = mpi::i_am_root() || !dist;
    size_t nmv_call = 0ul;
    bool stop = false;
    while (!stop){
        if (i_am_solver_rank) take_step();
        bool do_another = i_am_solver_rank && do_another_mv_call();
        if (dist) mpi::bcast(do_another);
        if (do_another){
            ++nmv_call;
            product_fn();
        }
        stop = i_am_solver_rank && basis_found();
        if (dist) mpi::bcast(stop);
    }
    // Finding eigenvalues and eigenvectors.
    if (i_am_solver_rank) find_eigenvalues();
    log::info("{} converged after {} {}distributed matrix-vector multiplication{}",
              string_utils::plural("eigenpair", m_nroot),
              nmv_call, dist ? "" : "non-",
              string_utils::plural(nmv_call));
    return true;
}
