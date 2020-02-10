//
// Created by Robert John Anderson on 2020-02-09.
//

#include <gtest/gtest.h>
#include <src/parallel/MPIWrapper.h>
#include <src/data/Specification.h>
#include "src/data/Table.h"

TEST(Table, AllToAllV) {
    MPIWrapper mpi;
    Specification spec;
    const size_t nint = 4;
    spec.create<size_t>(nint);
    const size_t nrow_send = 300;
    Table send_table(spec, nrow_send, mpi.nrank());
    Table recv_table(spec, nrow_send * mpi.nrank(), 1);

    // use the modulus to create a deterministic variability in the number of rows/entries
    auto scramble = [](const size_t isend, const size_t irecv, const size_t irow = 0,
                       const size_t ientry = 0, const size_t modular_divisor = ~0ul) {
        if (!modular_divisor) return 0ul;
        return ((isend + 9) * (isend * 9) +
                (irecv + 7) * (irecv + 7) +
                (irow + 5) * (irow + 5) +
                (ientry + 3) * (ientry + 3)) % modular_divisor;
    };

    for (auto isegment{0ul}; isegment < mpi.nrank(); ++isegment) {
        for (auto irow{0ul}; irow < scramble(mpi.irank(), isegment, 0, 0, nrow_send); ++irow) {
            send_table.claim_rows(isegment);
            for (auto iint{0ul}; iint < scramble(mpi.irank(), isegment, irow, 0, nint); ++iint) {
                *send_table.view<size_t>(isegment, irow, iint) =
                        scramble(mpi.irank(), isegment, irow, iint);
            }
        }
    }

    ASSERT_TRUE(send_table.send_to(recv_table));

    /*
     * check that all entries came through as expected
     */
    size_t irow_tot = 0;
    for (auto isend{0ul}; isend < mpi.nrank(); ++isend) {
        for (auto irow{0ul}; irow < scramble(isend, mpi.irank(), 0, 0, nrow_send); ++irow) {
            for (auto iint{0ul}; iint < scramble(isend, mpi.irank(), irow, 0, nint); ++iint) {
                ASSERT_EQ(
                        * recv_table.view<size_t>(0, irow_tot, iint),
                                scramble(isend, mpi.irank(), irow, iint)
                );
            }
            irow_tot++;
        }
    }
}