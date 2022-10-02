//
// Created by Robert J. Anderson on 10/11/2020.
//

#include <M7_lib/table/BufferedFields.h>
#include "M7_lib/communication/Communicator.h"
#include "gtest/gtest.h"

TEST(CommunicatingPair, CommunicateSingleElement) {
    typedef SingleFieldRow<field::Number<uint_t>> row_t;
    send_recv::BasicSend<row_t> comm_pair("Test pair", {}, {mpi::nrank(), 1.0});
    // after resize:
    const uint_t bw_size = comm_pair.row_size();
    const uint_t hash_lo = 123, hash_hi = 789;

    comm_pair.resize(1ul, 0.0);

    for (uint_t irank=0ul; irank<mpi::nrank(); ++irank) ASSERT_EQ(comm_pair.send(irank).nrecord(), 1ul);
    ASSERT_EQ(comm_pair.recv().nrecord(), mpi::nrank());

    for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto& row = comm_pair.send(irank).m_row;
        row.push_back_jump();
        row.m_field = hash::in_range({irank, mpi::irank()}, hash_lo, hash_hi);
    }

    for (uint_t idst = 1ul; idst < mpi::nrank(); ++idst) {
        // check pointer distance between adjacent send tables
        ASSERT_EQ(std::distance(comm_pair.send(idst - 1).begin(), comm_pair.send(idst).begin()), bw_size);
        // check send values by pointer dereference relative to entire send table array
        auto v_send = *reinterpret_cast<uint_t*>(comm_pair.send().begin() + idst*bw_size);
        ASSERT_EQ(v_send, hash::in_range({idst, mpi::irank()}, hash_lo, hash_hi));
    }

    comm_pair.communicate();

    auto& row = comm_pair.recv().m_row;
    for (row.restart(); row.in_range(); row.step()) {
        auto irank_src = row.index();
        ASSERT_EQ(row.m_field, hash::in_range({mpi::irank(), irank_src}, hash_lo, hash_hi));
    }
}

TEST(CommunicatingPair, CommunicateVectors){
    typedef SingleFieldRow<field::Numbers<uint_t, 1>> row_t;
    const double expansion_factor = 0.5;
    const uint_t nelement_vector = 13;
    send_recv::BasicSend<row_t> send_recv("Test pair", {nelement_vector}, {mpi::nrank(), expansion_factor});
    // after resize:
    const uint_t bw_size = send_recv.row_size();// * (1.0+expansion_factor);
    ASSERT_EQ(bw_size, nelement_vector*sizeof(uint_t));
    const uint_t hash_lo = 123, hash_hi = 789;

    send_recv.resize(1ul);
    for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto& row = send_recv.send(irank).m_row;
        row.push_back_jump();
        for (uint_t ielement=0ul; ielement<nelement_vector; ++ielement)
            row.m_field[ielement] = hash::in_range({irank, mpi::irank(), ielement}, hash_lo, hash_hi);
    }


    for (uint_t idst = 1ul; idst < mpi::nrank(); ++idst) {
        // check pointer distance between adjacent send tables
        ASSERT_EQ(std::distance(send_recv.send(idst - 1).begin(), send_recv.send(idst).begin()), bw_size);
        for (uint_t ielement=0ul; ielement<nelement_vector; ++ielement) {
            // check send values by pointer dereference relative to entire send table array
            auto v_send = reinterpret_cast<uint_t*>(send_recv.send().begin() + idst * bw_size)[ielement];
            ASSERT_EQ(v_send,hash::in_range({idst, mpi::irank(), ielement}, hash_lo, hash_hi));
        }
    }

    send_recv.communicate();
#if 0

    auto& row = send_recv.recv().m_row;
    for (row.restart(); row.in_range(); row.step()) {
        auto irank_src = row.index();
        for (uint_t ielement=0ul; ielement<nelement_vector; ++ielement)
            ASSERT_EQ(row.m_field[ielement],
                      hash::in_range({mpi::irank(), irank_src, ielement}, hash_lo, hash_hi));
    }
#endif
}

TEST(CommunicatingPair, CommunicateMultipleVectors){
    typedef SingleFieldRow<field::Numbers<uint_t, 1>> row_t;
    const double expansion_factor = 0.5;
    const uint_t nelement_vector = 5;
    const uint_t nrow_rank_lo = 6, nrow_rank_hi = 15;
    const uint_t nrow_this_rank = hash::in_range(mpi::irank(), nrow_rank_lo, nrow_rank_hi);
    send_recv::BasicSend<row_t> comm_pair("Test pair", {nelement_vector}, {nrow_this_rank, expansion_factor});
    const uint_t hash_lo = 123, hash_hi = 789;

    auto nrow_max = mpi::all_max(nrow_this_rank);

    v_t<uint_t> nrows_expect;
    nrows_expect.reserve(mpi::nrank());
    for (auto irank=0ul; irank<mpi::nrank(); ++irank)
        nrows_expect.push_back(hash::in_range(irank, nrow_rank_lo, nrow_rank_hi));
    ASSERT_EQ(*std::max_element(nrows_expect.cbegin(), nrows_expect.cend()), nrow_max);
    auto nrow_displs_expect = mpi::counts_to_displs_consec(nrows_expect);

    comm_pair.resize(nrow_max);
    for (uint_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        auto &row = comm_pair.send(irank).m_row;
        for (uint_t irow_send = 0ul; irow_send < nrow_this_rank; ++irow_send) {
            row.push_back_jump();
            for (uint_t ielement = 0ul; ielement < nelement_vector; ++ielement)
                row.m_field[ielement] = hash::in_range({irank, mpi::irank(), irow_send, ielement}, hash_lo, hash_hi);
        }
    }

    comm_pair.communicate();

    auto& row = comm_pair.recv().m_row;
    auto irank_src = 0ul;
    for (row.restart(); row.in_range(); row.step()) {
        if (irank_src+1<mpi::nrank()){
            if (row.index()==nrow_displs_expect[irank_src+1]) ++irank_src;
        }
        auto irow_send = row.index()-nrow_displs_expect[irank_src];
        for (uint_t ielement=0ul; ielement<nelement_vector; ++ielement)
            ASSERT_EQ(row.m_field[ielement],
                      hash::in_range({mpi::irank(), irank_src, irow_send, ielement}, hash_lo, hash_hi));
    }
}