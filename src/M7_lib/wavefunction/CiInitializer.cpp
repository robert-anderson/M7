//
// Created by anderson on 18/07/2022.
//

#include <M7_lib/util/ProgressMonitor.h>
#include "CiInitializer.h"
#include "M7_lib/foreach/ConnForeachGroup.h"
#include "M7_lib/field/Mbf.h"

ci_init::Subspace::Subspace(const Hamiltonian* h) :
        m_mbf_order_table("MBF order table", {mbf::row_t(h->m_basis, "mbf")}), m_h(h){
    m_mbf_order_table.set_expansion_factor(2.0);
}

ci_init::Initializer::Initializer(const Subspace& subspace, Options opts):
    m_opts(opts), m_is_hermitian(subspace.m_h->is_hermitian()){
    switch (opts.m_loop_kind) {
        case Options::Conns:
            build_ham_conns(subspace, opts.m_diag_shift);
            break;
        case Options::MbfPairs:
            build_ham_mbfs(subspace, opts.m_diag_shift);
            break;
    }
}

void ci_init::Initializer::build_ham_conns(const Subspace& subspace, ham_comp_t diag_shift) {
    const auto& table = subspace.m_mbf_order_table;
    const auto count = table.nrow_in_use();
    const uint_t count_local = mpi::evenly_shared_count(count);
    const uint_t displ_local = mpi::evenly_shared_displ(count);

    buffered::Mbf mbf(subspace.m_h->m_basis);
    conn::Mbf conn(mbf);
    ConnForeachGroup conn_iters(*subspace.m_h);

    m_sparse_ham.resize(count_local);

    logging::info("Building sparse H matrix ({} rows) by looping over MBFs and their connections", count_local);
    ProgressMonitor pm(true, "building sparse H", "basis functions", count_local, 5);
    auto& row = table.m_row;
    const auto& src_mbf = row.m_field;
    auto& dst_mbf = mbf;

    auto filling_fn = [&](){
        const auto helem = subspace.m_h->get_element(src_mbf, conn);
        if (!ham::is_significant(helem)) return;
        conn.apply(src_mbf, dst_mbf);
        auto& lookup = table.lookup(dst_mbf);
        DEBUG_ASSERT_TRUE(lookup, "connected MBF is outside generated space");
        const auto irow = row.index()-displ_local;
        if (lookup) m_sparse_ham.insert(irow, {lookup.index(), helem});
    };

    for (row.jump(displ_local); row.in_range(displ_local+count_local); ++row) {
        const auto irow = row.index() - displ_local;
        const auto helem_diag = subspace.m_h->get_element(src_mbf) + diag_shift;
        DEBUG_ASSERT_TRUE(m_sparse_ham[irow].empty(), "sparse row should be empty");
        if (ham::is_significant(helem_diag)) m_sparse_ham.insert(irow, {row.index(), helem_diag});
        src_mbf.m_decoded.clear();
        conn_iters.loop(conn, src_mbf, filling_fn);
        pm.next();
    }
}

void ci_init::Initializer::build_ham_mbfs(const Subspace& subspace, ham_comp_t diag_shift) {
    const auto& table = subspace.m_mbf_order_table;
    const auto count = table.nrow_in_use();
    const uint_t count_local = mpi::evenly_shared_count(count);
    const uint_t displ_local = mpi::evenly_shared_displ(count);

    m_sparse_ham.resize(count_local);
    conn::Mbf conn(subspace.m_h->m_basis.size());

    logging::info("Building sparse H matrix ({} rows) by looping over MBF pairs", count_local);
    ProgressMonitor pm(true, "building sparse H", "basis functions", count_local, 5);

    auto src = table.m_row;
    auto dst = src;

    for (src.jump(displ_local); src.in_range(displ_local + count_local); ++src) {
        const auto irow = src.index() - displ_local;
        DEBUG_ASSERT_TRUE(m_sparse_ham[irow].empty(), "sparse row should be empty");
        for (dst.restart(); dst; ++dst) {
            auto helem = subspace.m_h->get_element(src.m_field, dst.m_field);
            if (src.index() == dst.index()) helem += diag_shift;
            if (!ham::is_significant(helem)) continue;
            m_sparse_ham.insert(irow, {dst.index(), helem});
        }
        pm.next();
    }
}

ci_init::FciSubspace::FciSubspace(const Hamiltonian* h, sys::Particles particles) : Subspace(h) {
    auto iters = FciIters::make(*h, particles, false);
    m_mbf_order_table.resize(iters.niter_single());
    buffered::Mbf mbf(h->m_basis);
    iters.m_single->loop(mbf, [&](){m_mbf_order_table.insert(mbf);});
    m_mbf_order_table.remap();
}

ci_init::RefConnSubspace::RefConnSubspace(const Hamiltonian* h, const Mbf& ref) : Subspace(h) {
    conn::Mbf conn(ref);
    m_mbf_order_table.insert(ref);
    buffered::Mbf dst = ref;
    auto body = [&]() {
        DEBUG_ASSERT_NE(conn.exsig(), opsig::c_zero, "diagonal connection generated");
        auto helement = h->get_element(ref, conn);
        if (!ham::is_significant(helement)) return;
        conn.apply(ref, dst);
        m_mbf_order_table.insert(dst);
    };
    ConnForeachGroup(*h).loop(conn, ref, body);
    m_mbf_order_table.remap();
}