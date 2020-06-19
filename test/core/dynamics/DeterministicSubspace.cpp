//
// Created by rja on 08/06/2020.
//

#include <src/core/dynamics/WalkerList.h>
#include <src/core/hamiltonian/AbInitioHamiltonian.h>
#include <src/core/dynamics/DeterministicSubspace.h>
#include "gtest/gtest.h"

TEST(DeterministicSubspace, CreateSparseHamiltonian){
    ASSERT_TRUE(0);
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

    detsub.execute(1.0);
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