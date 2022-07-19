//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/util/ProgressMonitor.h>
#include "FciInitializer.h"

FciInitializer::FciInitializer(const Hamiltonian &h) {
    auto row_iters = FciIters::make(h);
    auto col_iters = FciIters::make(h);
    const auto count = row_iters.niter_single();

    auto count_local = mpi::evenly_shared_count(count);
    auto displ_local = mpi::evenly_shared_displ(count);

    buffered::Mbf row_mbf(h.m_basis);
    auto col_mbf = row_mbf;

    sparse::dynamic::Matrix<double> sparse_ham;
    sparse_ham.resize(count_local);

    ProgressMonitor pm(true, "building sparse H", "basis functions", count_local);
    uint_t irow = ~0ul;
    auto fn_row_loop = [&]() {
        ++irow;
        if ((irow < displ_local) || (irow >= (displ_local+count_local))) return;
        DEBUG_ASSERT_LT(irow-displ_local, uint_t(count_local), "local row index OOB");
        uint_t icol = ~0ul;
        auto fn_col_loop = [&]() {
            ++icol;
            auto helem = h.get_element(row_mbf, col_mbf);
            if (!ham::is_significant(helem)) return;
            sparse_ham.insert(irow-displ_local, {icol, helem});
        };
        col_iters.m_single->loop(col_mbf, fn_col_loop);
        pm.next();
    };
    row_iters.m_single->loop(row_mbf, fn_row_loop);

    const uint_t nroot = 3;
    ArnoldiProblemNonSym<double> solver(nroot);
    dist_mv_prod::Sparse<double> dist(sparse_ham);
    solver.solve(dist);
    if (mpi::i_am_root()) m_eval = solver.real_eigenvalue(nroot-1);
    mpi::bcast(m_eval);
    logging::info("Non-symmetric Arnoldi eigenvalue: {}", m_eval);
}
