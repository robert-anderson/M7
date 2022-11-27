//
// Created by Robert J. Anderson on 09/01/2022.
//

#include "ArnoldiSolver.h"
#include "M7_lib/util/String.h"


bool ArnoldiSolverBase::solve(const std::function<void()> &product_fn, bool dist) {
    bool i_am_solver_rank = mpi::i_am_root() || !dist;
    uint_t nmv_call = 0ul;
    char stop = false;
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
    char success = false;
    if (i_am_solver_rank) {
        success = find_eigenvectors();
        logging::info("ARPACK called {}distributed matrix-vector multiplication {} time{}",
                      (dist ? "" : "non-"), nmv_call, string::plural(nmv_call));
    }
    if (dist) mpi::bcast(success);
    return success;
}