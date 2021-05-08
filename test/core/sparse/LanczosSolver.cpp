//
// Created by rja on 07/05/2021.
//

#include <src/core/table/BufferedFields.h>
#include <src/core/hamiltonian/FermionHamiltonian.h>
#include "gtest/gtest.h"
#include "src/core/sparse/LanczosSolver.h"


struct ProgressMonitor {
    const bool m_local;
    const std::string m_name, m_item_name;
    const size_t m_nexpect, m_pc_resolution, m_period;
    size_t m_i = 0ul;

    ProgressMonitor(bool local, std::string name, std::string item_name, size_t nexpect, size_t pc_resolution = 5) :
            m_local(local), m_name(name), m_item_name(item_name), m_nexpect(nexpect),
            m_pc_resolution(pc_resolution), m_period(std::round(nexpect * (pc_resolution / 100.0))) {
        if (m_local)
            log::info_("Starting process \"{}\"...", m_name);
        else
            log::info("Starting process \"{}\"...", m_name);
    }

private:
    void log(const size_t &pc) const {
        if (m_local)
            log::info_("process \"{}\" is {}% ({}/{} {}) complete",
                       m_name, pc, m_i + 1, m_nexpect, m_item_name);
        else
            log::info("process \"{}\" is {}% ({}/{} {}) complete",
                   m_name, pc, m_i + 1, m_nexpect, m_item_name);
    }

public:
    void next() {
        auto nperiod = m_i / m_period;
        if (m_i + 1 == m_nexpect) log(100);
        else if (m_i && !(m_i % m_period)) log(nperiod * m_pc_resolution);
        ++m_i;
    }
};

TEST(LanczosSolver, Test) {
    typedef SingletRow<fields::FermionOnv> onv_row_t;
    FermionHamiltonian ham(defs::assets_root + "/Hubbard_U4_8site/FCIDUMP", false);
    BufferedTable<onv_row_t, true> onv_table("", {{ham.nsite()}, 800000});
    buffered::FermionOnv work_onv(ham.nsite());
    enums::CombinationsDistinct alpha_combs(ham.nsite(), ham.nelec() / 2);

    auto nci_1spin = integer_utils::combinatorial(ham.nsite(), ham.nelec() / 2);
    ProgressMonitor build_fci_progress(true, "build FCI space", "alpha combinations", nci_1spin);
    while (alpha_combs.next()) {
        build_fci_progress.next();
        enums::CombinationsDistinct beta_combs(ham.nsite(), ham.nelec() / 2);
        while (beta_combs.next()) {
            work_onv.set(alpha_combs.m_v, beta_combs.m_v);
            onv_table.insert(work_onv);
        }
    }
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