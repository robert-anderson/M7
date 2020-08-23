//
// Created by rja on 08/06/2020.
//

#include <src/core/dynamics/WalkerList.h>
#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include <src/core/dynamics/DeterministicSubspace.h>
#include <src/core/util/Timer.h>
#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include "gtest/gtest.h"

TEST(DeterministicSubspace, FciCheck) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    auto ref_det = ham.guess_reference(0);
    WalkerList walker_list("test walker list", ham.nsite(), 100);
    RankAllocator<DeterminantElement> ra(10, 1);
    ham.generate_ci_space(&walker_list, ra, 0);
    ASSERT_EQ(mpi::all_sum(walker_list.high_water_mark(0)), 400);
    bool is_ref_rank = walker_list.m_determinant(0) == ref_det;
    if (is_ref_rank)
        walker_list.m_weight(0) = 1;
    DeterministicSubspace detsub(walker_list);
    detsub.build_from_whole_walker_list(&ham);
    double tau = 0.02;

    auto do_iter = [&]() {
        Reducible<defs::ham_t> num;
        Reducible<defs::ham_comp_t> norm_square;
        Reducible<defs::wf_comp_t> delta_nw;
        detsub.gather_and_project();
        detsub.rayleigh_quotient(num, norm_square);
        auto l1_init = walker_list.l1_norm(0);
        for (size_t i = 0; i < walker_list.high_water_mark(0); ++i) {
            auto hdiag = *walker_list.m_hdiag(i);
            auto weight = *walker_list.m_weight(i);
            num.local() += hdiag * std::pow(std::abs(weight), 2.0);
            delta_nw.local() -= std::abs(*walker_list.m_weight(i));
            walker_list.m_weight(i) *= 1 - tau * hdiag;
            delta_nw.local() += std::abs(*walker_list.m_weight(i));
        }
        detsub.update_weights(tau, delta_nw);
        auto l1_final = walker_list.l1_norm(0);
        delta_nw.mpi_sum();
        ASSERT(consts::floats_nearly_equal(l1_init + delta_nw.reduced(), l1_final));
        walker_list.normalize();
        return consts::real(num.mpi_sum()) / norm_square.mpi_sum();
    };

    defs::ham_comp_t energy;
    energy = do_iter();
    if (is_ref_rank) {
        ASSERT_FLOAT_EQ(energy, ham.get_energy(ref_det));
    }

    energy = do_iter();
    ASSERT_FLOAT_EQ(energy, -108.65403);

    const size_t niter = 10000;
    for (size_t i = 0ul; i < niter; ++i) {
        energy = do_iter();
    }
    ASSERT_TRUE(consts::floats_nearly_equal(energy, -108.8113865756313, 1e-6));
}

TEST(DeterministicSubspace, BuildFromDeterminantConnections) {
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    WalkerList walker_list("test walker list", ham.nsite(), 100);
    auto ref = ham.guess_reference(0);
    auto excited = ref;
    const size_t nconn = 48;

    HamiltonianConnectionEnumerator enumerator(ham, ref);
    MatrixElement<defs::ham_t> matel(ref);

    ASSERT_LT(mpi::nrank(), nconn);
    const size_t nconn_per_rank = nconn / mpi::nrank();
    walker_list.expand(nconn_per_rank);
    for (size_t iconn = 0ul; iconn < nconn; ++iconn) {
        // add to walker list on this process
        enumerator.next(matel);
        if (iconn%mpi::nrank()!=mpi::irank()) continue;
        matel.aconn.apply(ref, excited);
        auto irow = walker_list.push(excited);
        walker_list.m_weight(irow) = mpi::irank() + 1;
    }
    DeterministicSubspace detsub(walker_list);
    detsub.build_from_det_connections(ref, &ham);

    ASSERT_EQ(detsub.nrow_local(), walker_list.high_water_mark(0));
    ASSERT_EQ(detsub.nrow_full(), mpi::nrank() * nconn_per_rank);

    detsub.gather_and_project();
    Reducible<defs::wf_comp_t> delta_nw;
    detsub.update_weights(1.0, delta_nw);
}

TEST(DeterministicSubspace, BuildFromHighestWeighted) {
    /*
    AbInitioHamiltonian ham(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_TRUE(ham.spin_conserving());
    WalkerList walker_list("test walker list", ham.nsite(), 100);
    auto ref = ham.guess_reference(0);
    DeterministicSubspace detsub(walker_list);
    detsub.build_from_det_connections(ref, &ham);

    ASSERT_EQ(detsub.nrow_local(), walker_list.high_water_mark(0));
    ASSERT_EQ(detsub.nrow_full(), mpi::nrank() * nconn_per_rank);

    detsub.gather_and_project();
    Reducible<defs::wf_comp_t> delta_nw;
    detsub.update_weights(1.0, delta_nw);
     */
}

TEST(DeterministicSubspace, SerializeToDisk) {
    ASSERT_TRUE(0);
}

TEST(DeterministicSubspace, BuildFromSerialized) {
    ASSERT_TRUE(0);
}
