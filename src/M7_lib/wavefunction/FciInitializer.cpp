//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/util/ProgressMonitor.h>
#include "FciInitializer.h"
#include "M7_lib/foreach/ConnForeachGroup.h"

FciInitializer::FciInitializer(const Hamiltonian &h, FciInitOptions opts):
        m_mbf_order_table("MBF order table", {{h.m_basis}}){
    auto iters = FciIters::make(h);
    const auto count = iters.niter_single();
    const uint_t count_local = mpi::evenly_shared_count(count);
    const uint_t displ_local = mpi::evenly_shared_displ(count);

    m_mbf_order_table.resize(iters.niter_single());
    buffered::Mbf mbf(h.m_basis);
    auto setup_fn = [&](){m_mbf_order_table.insert(mbf);};
    iters.m_single->loop(mbf, setup_fn);

    conn::Mbf conn(mbf);
    ConnForeachGroup conn_iters(h);

    sparse::dynamic::Matrix<double> sparse_ham;
    sparse_ham.resize(count_local);

    logging::info("Building sparse H matrix ({} rows)", iters.niter_single());
    ProgressMonitor pm(true, "building sparse H", "basis functions", count_local);
    auto& row = m_mbf_order_table.m_row;
    for (row.jump(displ_local); row.in_range(displ_local+count_local); row.step()) {
        const auto& src_mbf = row.m_mbf;
        auto& dst_mbf = mbf;
        const auto irow = row.index()-displ_local;

        const auto helem_diag = h.get_element(src_mbf) + opts.m_diag_shift;
        DEBUG_ASSERT_TRUE(sparse_ham[irow].empty(), "sparse row should be empty");
        if (ham::is_significant(helem_diag)) sparse_ham.insert(irow, {row.index(), helem_diag});
        auto filling_fn = [&](){
            const auto helem = h.get_element(src_mbf, conn);
            if (!ham::is_significant(helem)) return;
            conn.apply(src_mbf, dst_mbf);
            const auto icol = m_mbf_order_table.lookup(dst_mbf).index();
            sparse_ham.insert(irow, {icol, helem});
        };
        src_mbf.m_decoded.clear();
        conn_iters.loop(conn, src_mbf, filling_fn);
        pm.next();
    }

    /*
     * once the ARPACK procedure is complete, the eigenvalues must be adjusted to undo the diagonal shift
     */
    auto undo_shift = [&opts](ARrcStdEig<ham_t, ham_t>* arpack) {
        if (!arpack) return;
        for (auto ptr = arpack->RawEigenvalues(); ptr!=arpack->RawEigenvalues()+arpack->GetNev(); ++ptr)
            *ptr-=opts.m_diag_shift;
    };

    dist_mv_prod::Sparse<double> dist(sparse_ham);
    if (h.is_hermitian()) {
        m_arpack_sym.solve(dist, opts);
        undo_shift(m_arpack_sym.m_solver.get());
    }
    else {
        m_arpack_nonsym.solve(dist, opts);
        undo_shift(m_arpack_nonsym.m_solver.get());
    }
}
