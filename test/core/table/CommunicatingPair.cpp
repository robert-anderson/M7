//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/table/Communicator.h"

#if 0
TEST(CommunicatingPair, CommunicateSingleElement) {
    struct TestTable : public Table {
        fields::Number<size_t> m_counter;

        TestTable() : m_counter(this, "counter") {}
    };

    double expansion_factor = 0.5;
    CommunicatingPair<TestTable> comm_pair("Test pair", expansion_factor, {});
    // after resize:
    const size_t ndword_bw = comm_pair.row_dsize()*mpi::nrank();

    comm_pair.resize(mpi::nrank());
    for (size_t irank = 0ul; irank < mpi::nrank(); ++irank) {
        comm_pair.send(irank).push_back();
        comm_pair.send(irank).m_counter(0) = 9900 + 10 * irank + mpi::irank();
    }

    for (size_t idst = 1ul; idst < mpi::nrank(); ++idst) {
        ASSERT_EQ(std::distance(comm_pair.send(idst - 1).dbegin(), comm_pair.send(idst).dbegin()), ndword_bw);
        ASSERT_EQ(comm_pair.send().dbegin()[idst*ndword_bw], 9900 + 10 * idst + mpi::irank());
    }

    comm_pair.communicate();

    for (size_t isrc = 0ul; isrc < mpi::nrank(); ++isrc)
        ASSERT_EQ(comm_pair.recv().m_counter(isrc), 9900 + 10 * mpi::irank() + isrc);
}

TEST(CommunicatingPair, CommunicateVector) {
    struct TestTable : public Table {
        fields::Numbers<int, 1> m_counter;
        TestTable(size_t nint) : m_counter(this, "counter", nint) {}
    };

    const size_t nrow = 120;
    const size_t nint = 7;

    CommunicatingPair<TestTable> comm_pair("Test pair", 0.5, {nint});

    auto flat_index = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                         const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return (((isend * nrow + irecv) * nrow + irow * nrow) + ientry) % modular_divisor;
    };

    comm_pair.resize(nrow);

    //ASSERT_EQ(comm_pair.recv().bw_dsize(), mpi::nrank() * nrow_ * comm_pair.row_dsize());
    ASSERT_EQ(comm_pair.send().buffer_dsize(), mpi::nrank() * nrow * comm_pair.row_dsize());

    for (size_t idst_rank = 0ul; idst_rank < mpi::nrank(); ++idst_rank) {
        for (size_t irow = 0ul; irow < nrow; ++irow) {
            ASSERT_EQ(comm_pair.send(idst_rank).push_back(), irow);
            for (size_t iint = 0ul; iint < nint; ++iint) {
                comm_pair.send(idst_rank).m_counter(irow, iint) = flat_index(mpi::irank(), idst_rank, irow, iint);
            }
        }
    }

    comm_pair.communicate();

/*
 * check that all entries came through as expected
 */
    size_t irow_tot = 0;
    for (size_t isend = 0ul; isend < mpi::nrank(); ++isend) {
        for (size_t irow = 0ul; irow < nrow; ++irow) {
            for (size_t iint = 0ul; iint < nint; ++iint) {
                ASSERT_EQ(
                        comm_pair.recv().m_counter(irow_tot, iint),
                        flat_index(isend, mpi::irank(), irow, iint)
                );
            }
            irow_tot++;
        }
    }
}
#endif