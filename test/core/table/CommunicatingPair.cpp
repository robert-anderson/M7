//
// Created by rja on 10/11/2020.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/table/Communicator.h"
#include "gtest/gtest.h"

TEST(CommunicatingPair, CommunicateSingleElement) {
    typedef SingleFieldRow<fields::Number<size_t>> row_t;
    const double expansion_factor = 0.5;
    CommunicatingPair<row_t> comm_pair("Test pair", expansion_factor, {{}});
    // after resize:
    const size_t ndword_bw = comm_pair.row_dsize()*mpi::nrank();
    const size_t hash_lo = 123, hash_hi = 789;

    comm_pair.resize(mpi::nrank());
    for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto& row = comm_pair.send(irank).m_row;
        row.push_back_jump();
        row.m_field = hashing::in_range({irank, mpi::irank()}, hash_lo, hash_hi);
    }

    for (size_t idst = 1ul; idst < mpi::nrank(); ++idst) {
        // check pointer distance between adjacent send tables
        ASSERT_EQ(std::distance(comm_pair.send(idst - 1).dbegin(), comm_pair.send(idst).dbegin()), ndword_bw);
        // check send values by pointer dereference relative to entire send table array
        ASSERT_EQ(comm_pair.send().dbegin()[idst*ndword_bw], hashing::in_range({idst, mpi::irank()}, hash_lo, hash_hi));
    }

    comm_pair.communicate();

    auto& row = comm_pair.recv().m_row;
    for (row.restart(); row.in_range(); row.step()) {
        auto irank_src = row.m_i;
        ASSERT_EQ(row.m_field, hashing::in_range({mpi::irank(), irank_src}, hash_lo, hash_hi));
    }
}

TEST(CommunicatingPair, CommunicateVectors){
    typedef SingleFieldRow<fields::Numbers<size_t, 1>> row_t;
    const double expansion_factor = 0.5;
    const size_t nelement_vector = 13;
    CommunicatingPair<row_t> comm_pair("Test pair", expansion_factor, {{nelement_vector}});
    // after resize:
    const size_t ndword_bw = comm_pair.row_dsize()*mpi::nrank();
    const size_t hash_lo = 123, hash_hi = 789;

    comm_pair.resize(mpi::nrank());
    for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto& row = comm_pair.send(irank).m_row;
        row.push_back_jump();
        for (size_t ielement=0ul; ielement<nelement_vector; ++ielement)
            row.m_field[ielement] = hashing::in_range({irank, mpi::irank(), ielement}, hash_lo, hash_hi);
    }

    for (size_t idst = 1ul; idst < mpi::nrank(); ++idst) {
        // check pointer distance between adjacent send tables
        ASSERT_EQ(std::distance(comm_pair.send(idst - 1).dbegin(), comm_pair.send(idst).dbegin()), ndword_bw);
        for (size_t ielement=0ul; ielement<nelement_vector; ++ielement)
            // check send values by pointer dereference relative to entire send table array
            ASSERT_EQ(comm_pair.send().dbegin()[idst*ndword_bw+ielement],
                      hashing::in_range({idst, mpi::irank(), ielement}, hash_lo, hash_hi));
    }

    comm_pair.communicate();

    auto& row = comm_pair.recv().m_row;
    for (row.restart(); row.in_range(); row.step()) {
        auto irank_src = row.m_i;
        for (size_t ielement=0ul; ielement<nelement_vector; ++ielement)
            ASSERT_EQ(row.m_field[ielement],
                      hashing::in_range({mpi::irank(), irank_src, ielement}, hash_lo, hash_hi));
    }
}

TEST(CommunicatingPair, CommunicateMultipleVectors){
    typedef SingleFieldRow<fields::Numbers<size_t, 1>> row_t;
    const double expansion_factor = 0.5;
    const size_t nelement_vector = 5;
    const size_t nrow_rank_lo = 6, nrow_rank_hi = 15;
    const size_t nrow_this_rank = hashing::in_range(mpi::irank(), nrow_rank_lo, nrow_rank_hi);
    CommunicatingPair<row_t> comm_pair("Test pair", expansion_factor, {{nelement_vector}});
    const size_t hash_lo = 123, hash_hi = 789;

    auto nrow_max = mpi::all_max(nrow_this_rank);

    std::vector<size_t> nrows_expect;
    nrows_expect.reserve(mpi::nrank());
    for (auto irank=0ul; irank<mpi::nrank(); ++irank)
        nrows_expect.push_back(hashing::in_range(irank, nrow_rank_lo, nrow_rank_hi));
    ASSERT_EQ(*std::max_element(nrows_expect.cbegin(), nrows_expect.cend()), nrow_max);
    auto nrow_displs_expect = mpi::counts_to_displs_consec(nrows_expect);

    comm_pair.resize(nrow_max);
    for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto &row = comm_pair.send(irank).m_row;
        for (size_t irow_send = 0ul; irow_send < nrow_this_rank; ++irow_send) {
            row.push_back_jump();
            for (size_t ielement = 0ul; ielement < nelement_vector; ++ielement)
                row.m_field[ielement] = hashing::in_range({irank, mpi::irank(), irow_send, ielement}, hash_lo, hash_hi);
        }
    }

    comm_pair.communicate();

    auto& row = comm_pair.recv().m_row;
    auto irank_src = 0ul;
    for (row.restart(); row.in_range(); row.step()) {
        if (irank_src+1<mpi::nrank()){
            if (row.m_i==nrow_displs_expect[irank_src+1]) ++irank_src;
        }
        auto irow_send = row.m_i-nrow_displs_expect[irank_src];
        for (size_t ielement=0ul; ielement<nelement_vector; ++ielement)
            ASSERT_EQ(row.m_field[ielement],
                      hashing::in_range({mpi::irank(), irank_src, irow_send, ielement}, hash_lo, hash_hi));
    }
}