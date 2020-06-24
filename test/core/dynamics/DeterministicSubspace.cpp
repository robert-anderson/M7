//
// Created by rja on 08/06/2020.
//

#include <src/core/dynamics/WalkerList.h>
#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include <src/core/dynamics/DeterministicSubspace.h>
#include "gtest/gtest.h"

TEST(DeterministicSubspace, FciCheck){
    AbInitioHamiltonian ham(defs::assets_root+"/RHF_N2_6o6e/FCIDUMP");
    auto ref_det = ham.guess_reference(0);
    WalkerList walker_list(ham.nsite(), 100);
    RankAllocator<DeterminantElement> ra(10);
    ham.generate_ci_space(&walker_list, ra, 0);
    ASSERT_EQ(mpi::all_sum(walker_list.high_water_mark(0)), 400);
    bool is_ref_rank = walker_list.m_determinant(0)==ref_det;
    if (is_ref_rank)
        walker_list.m_weight(0)=1;
    DeterministicSubspace detsub(&walker_list);
    detsub.build_from_whole_walker_list(&ham);
    double tau = 0.02;

    auto do_iter = [&](){
        Hybrid<defs::ham_t> num;
        Hybrid<defs::ham_comp_t> norm_square;
        detsub.gather_and_project();
        detsub.rayleigh_quotient(num, norm_square);
        for (size_t i=0; i< walker_list.high_water_mark(0); ++i){
            auto hdiag = *walker_list.m_hdiag(i);
            auto weight = *walker_list.m_weight(i);
            num.local()+=hdiag*std::pow(std::abs(weight), 2.0);
            walker_list.m_weight(i)*=1-tau*hdiag;
        }
        detsub.update_weights(tau);
        walker_list.normalize();
        return consts::real(num.mpi_sum())/norm_square.mpi_sum();
    };

    defs::ham_comp_t energy;
    energy = do_iter();
    if (is_ref_rank)
        ASSERT_FLOAT_EQ(energy, ham.get_energy(ref_det));

    energy = do_iter();
    ASSERT_FLOAT_EQ(energy, -108.65403);

    const size_t niter = 10000;
    for (size_t i=0ul; i<niter; ++i){
        energy = do_iter();
    }
    ASSERT_TRUE(consts::floats_nearly_equal(energy, -108.8113865756313, 1e-6));
}

TEST(DeterministicSubspace, BuildFromDeterminantConnections){
    AbInitioHamiltonian ham(defs::assets_root+"/RHF_N2_6o6e/FCIDUMP");
    WalkerList walker_list(ham.nsite(), 100);
    auto ref = ham.guess_reference(0);
    Hamiltonian::ConnectionList conn_list(ham.nsite(), 100);
    const size_t nconn = 48;
    conn_list.expand(nconn);
    auto irow_ref = conn_list.push(ref);

    ham.all_connections_of_det(&conn_list, ref, 1e-12);
    ASSERT_EQ(conn_list.high_water_mark(0), nconn);
    ASSERT_LT(mpi::nrank(), nconn);
    const size_t nconn_per_rank = nconn/mpi::nrank();
    walker_list.expand(nconn_per_rank);
    for (size_t iconn=mpi::irank()*nconn_per_rank; iconn<(mpi::irank()+1)*nconn_per_rank; ++iconn){
        // add to walker list on this process
        auto det = conn_list.determinant(iconn);
        auto irow = walker_list.push(det);
        walker_list.m_weight(irow) = mpi::irank()+1;
    }
    DeterministicSubspace detsub(&walker_list);
    detsub.build_from_det_connections(ref, &ham);

    ASSERT_EQ(detsub.nrow_local(), walker_list.high_water_mark(0));
    ASSERT_EQ(detsub.nrow_full(), mpi::nrank()*nconn_per_rank);

    detsub.gather_and_project();
    detsub.update_weights(1.0);

}

TEST(DeterministicSubspace, BuildFromHighestWeighted){
    ASSERT_TRUE(0);
}

TEST(DeterministicSubspace, SerializeToDisk){
    ASSERT_TRUE(0);
}

TEST(DeterministicSubspace, BuildFromSerialized){
    ASSERT_TRUE(0);
}