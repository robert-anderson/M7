//
// Created by Robert John Anderson on 2020-02-09.
//

#include <gtest/gtest.h>
#include "src/data/Table.h"
#include "src/enumerators/BitfieldEnumerator.h"

TEST(Table, EncodeDecode) {
    struct TestTable : public Table {
        Field<BitfieldNew> m_bitfield;
        Field<size_t> m_integers;
        TestTable(size_t nbit, size_t nint):
        Table(),
        m_bitfield(this, nbit),
        m_integers(this, nint){}
    };
    const size_t nint = 4;
    const size_t nbit = 20;

    TestTable table(nbit, nint);
    table.grow(8);
    table.m_integers.get(3, 3) = 6;
    defs::inds setinds = {3, 7, 14, 16, 19};
    table.m_bitfield.get(7).set(setinds);
    table.print();
    ASSERT_EQ(table.m_integers.get(3, 3), 6);
    auto setind = setinds.begin();
    BitfieldSetEnumerator enumerator(table.m_bitfield.get(7));
    size_t encoded_setind;
    while (enumerator.next(encoded_setind)){
        ASSERT_EQ(encoded_setind, *(setind++));
    }
}

/*
TEST(Table, AllToAllV) {
    Specification spec;
    const size_t nint = 4;
    spec.add<size_t>(nint);
    const size_t nrow_send = 10;
    Table send_table(spec, nrow_send, mpi::nrank());
    Table recv_table(spec, nrow_send * mpi::nrank(), 1);

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
    size_t irow_tot = 0;
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
}*/

    /*

TEST(Table, PrivateTempThreadSafety) {
    Specification specification;
    specification.add<size_t>(3);
    size_t nrow = 3600;
    Table shared_table(specification, nrow);

    auto f = [](size_t i, size_t ithread) {
        return integer_utils::combinatorial(2 + ithread + i % 7, 2);
    };

#pragma omp parallel
    {
        Table private_table(shared_table, 7);
#pragma omp for
        for (size_t i = 0; i < nrow; ++i) {
            size_t irow = private_table.push(0, 1);
            private_table.view<size_t>(irow)[0] = i;
            private_table.view<size_t>(irow)[1] = omp_get_thread_num();
            private_table.view<size_t>(irow)[2] = f(i, omp_get_thread_num());
        }
    }

    for (auto irow{0ul}; irow < shared_table.nrow(); ++irow) {
        auto v = shared_table.view<size_t>(irow);
        ASSERT_EQ(v[2], f(v[0], v[1]));
    }
}*/