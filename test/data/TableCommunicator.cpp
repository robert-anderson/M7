//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#include <src/data/Specification.h>
#include "src/data/Table.h"
#include "src/data/TableCommunicator.h"
#include "src/data/BitfieldNew.h"

/*
TEST(TableCommunicator, AllToAllV) {
    Specification spec;
    const size_t nint = 4;
    spec.add<size_t>(nint);
    const size_t nrow_send = 10;
    const size_t nrow_recv = nrow_send*mpi::nrank();

    TableCommunicator<Table> communicator(spec, nrow_send, nrow_recv);

    auto flat_index = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                         const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return (((isend * nrow_send + irecv) * nrow_send + irow * nrow_send) + ientry) % modular_divisor;
    };

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
    */
    /*
     * check that all entries came through as expected
     */
    /*
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
};
     */