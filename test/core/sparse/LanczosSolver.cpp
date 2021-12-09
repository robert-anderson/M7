//
// Created by rja on 07/05/2021.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/hamiltonian/FrmHam.h>
#include <src/core/basis/CiSpaces.h>
#include <src/core/util/ProgressMonitor.h>
#include "gtest/gtest.h"

#if 0
#include "src/core/sparse/LanczosSolver.h"


#ifndef ENABLE_BOSONS
TEST(LanczosSolver, Test) {
    typedef SingleFieldRow<fields::FermionOnv> onv_row_t;
    FermionHamiltonian ham(defs::assets_root + "/Hubbard_U4_8site/FCIDUMP", false);
    BufferedTable<onv_row_t, true> onv_table("", {{ham.nsite()}, 800000});
    enums::CombinationsDistinct alpha_combs(ham.nsite(), ham.nelec() / 2);

    auto nci_1spin = integer_utils::combinatorial(ham.nsite(), ham.nelec() / 2);
    //ProgressMonitor build_fci_progress(true, "build FCI space", "alpha combinations", nci_1spin);
    ci_gen::SpinSym gen(ham.nsite(), ham.nelec(), 0);
    gen(onv_table.m_row, onv_table.m_row.m_field);
    ASSERT_EQ(onv_table.m_hwm, nci_1spin * nci_1spin);

    conn::Antisym<0> conn(ham.nsite());
    sparse::Matrix<defs::ham_t> sparse_mat;
    sparse_mat.resize(onv_table.m_hwm);
    auto row = onv_table.m_row;
    auto row2 = onv_table.m_row;

    auto fn_body = [&](const conn::Antisym<0> &conn, const fields::FermionOnv &dst_onv, const defs::ham_t &helem) {
        ASSERT(!consts::float_is_zero(helem));
        auto irow_dst = *onv_table[dst_onv];
        sparse_mat.add(row.m_i, irow_dst, helem);
    };

    ProgressMonitor build_ham_progress(true, "build sparse Hamiltonian", "rows", onv_table.m_hwm, 3);
    for (row.restart(); row.in_range(); row.step()) {
        build_ham_progress.next();
        ham.foreach_connection(row.m_field, fn_body, true, true, true);
    }

    LanczosSolver solver(1);
    solver.solve(sparse_mat, 40);
    std::cout << solver.m_evals[0] << std::endl;
}
#endif
#endif