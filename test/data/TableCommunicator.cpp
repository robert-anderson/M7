//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#if 0
#include <src/data/Specification.h>
#include "src/data/Table.h"
#include "src/data/TableCommunicator.h"
#include "src/data/BitfieldNew.h"

TEST(TableCommunicator, AllToAllV) {
    struct TestTable : public Table {
        Field<size_t> m_integers;
        TestTable(size_t nint): Table(),
        m_integers(this, nint){}
    };
    const size_t nint = 4;
    TestTable table(nint);

    const size_t nrow = 6;

    auto factory = [nint](){return TestTable(nint);};
    TableCommunicator<TestTable> communicator(factory, nrow);

    auto flat_index = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                         const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return (((isend * nrow + irecv) * nrow + irow * nrow) + ientry) % modular_divisor;
    };

    for (auto idst_rank{0ul}; idst_rank < mpi::nrank(); ++idst_rank) {
        for (auto irow{0ul}; irow < nrow; ++irow) {
            auto send_table = communicator.m_send[idst_rank];
            for (auto iint{0ul}; iint < nint; ++iint) {
                send_table.m_integers.get(irow, iint) = flat_index(mpi::irank(), idst_rank, irow, iint);
            }
        }
    }

    communicator.m_send.print();


#if 0


    for (auto isegment{0ul}; isegment < mpi::nrank(); ++isegment) {
        for (auto irow{0ul}; irow < nrow_send; ++irow) {
            send_table.push(isegment, 1);
            for (auto iint{0ul}; iint < nint; ++iint) {
                *send_table.view<size_t>(isegment, irow, iint) =
                    flat_index(mpi::irank(), isegment, irow, iint);
            }
        }
    }
    ASSERT_TRUE(send_table.send_to(recv_table));
    /*
     * check that all entries came through as expected
     */
    size_t
    irow_tot = 0;
    for (auto isend{0ul}; isend < mpi::nrank(); ++isend) {
        for (auto irow{0ul}; irow < nrow_send; ++irow) {
            for (auto iint{0ul}; iint < nint; ++iint) {
                ASSERT_EQ(
                    *recv_table.view<size_t>(0, irow_tot, iint),
                    flat_index(isend, mpi::irank(), irow, iint)
                );
            }
            irow_tot++;
        }
    }
#endif
}
#endif