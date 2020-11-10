//
// Created by rja on 10/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/table/Table.h"
#include "src/core/table/CommunicatingPair.h"

struct TestTable : public TableX {
    fields::Numbers<int, 1> m_counter;
    TestTable(size_t nint): m_counter(this, "counter", nint){}
};


TEST(CommunicatingPair, Communication) {
    const size_t nrow = 100;
    const size_t nint = 6;

    CommunicatingPair<TestTable> comm_pair(nint);

    auto flat_index = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                         const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return (((isend * nrow + irecv) * nrow + irow * nrow) + ientry) % modular_divisor;
    };

    comm_pair.expand(nrow);

    ASSERT(comm_pair.recv().bw_dsize() == nrow*comm_pair.row_dsize())
    ASSERT(comm_pair.send().buffer_dsize() == nrow*comm_pair.row_dsize())

    for (auto idst_rank{0ul}; idst_rank < mpi::nrank(); ++idst_rank) {
        for (auto irow{0ul}; irow < nrow; ++irow) {
            ASSERT_EQ(comm_pair.send(idst_rank).push_back(), irow);
            for (auto iint{0ul}; iint < nint; ++iint) {
                comm_pair.send(idst_rank).m_counter(irow, iint) = flat_index(mpi::irank(), idst_rank, irow, iint);
            }
        }
    }

    comm_pair.communicate();
/*
 * check that all entries came through as expected
 */
    size_t irow_tot = 0;
    for (auto isend{0ul}; isend < mpi::nrank(); ++isend) {
        for (auto irow{0ul}; irow < nrow; ++irow) {
            for (auto iint{0ul}; iint < nint; ++iint) {
                ASSERT_EQ(
                        comm_pair.recv().m_counter(irow_tot, iint),
                        flat_index(isend, mpi::irank(), irow, iint)
                );
            }
            irow_tot++;
        }
    }
}
