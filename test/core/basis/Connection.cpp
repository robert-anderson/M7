//
// Created by Robert John Anderson on 2020-03-31.
//

#include "gtest/gtest.h"
#include "src/core/basis/FermionOnvConnection.h"
#include "src/core/io/SparseArrayFileReader.h"
#include "src/core/table/BufferedFields.h"
#include "src/core/hamiltonian/Hamiltonian.h"


TEST(Connection, ParticleNumberConserving) {
    const size_t nsite = 70;
    buffered::FermionOnv ket(nsite);
    buffered::FermionOnv bra(nsite);

    defs::inds ketoccorbs = {1, 4, 6, 8, 11, 19, 120, 138, 139};
    defs::inds braoccorbs = {1, 4, 5, 6, 9, 11, 19, 137, 138};
    ASSERT_EQ(ket.m_size, 3 * defs::nbyte_word);
    ASSERT_EQ(ket.m_dsize, 3);
    ASSERT_EQ(ket.m_nbit_in_last_dword, nsite * 2 - 2 * 64);
    ket = ketoccorbs;
    bra = braoccorbs;

    ASSERT_EQ(ket.nsetbit(), 9);
    ASSERT_EQ(bra.nsetbit(), 9);

    FermionOnvConnection conn(ket, bra);

    ASSERT_EQ(conn.nann(), 3);
    ASSERT_EQ(conn.ncre(), 3);

    ASSERT_EQ(conn.ann(0), 8);
    ASSERT_EQ(conn.ann(1), 120);
    ASSERT_EQ(conn.ann(2), 139);

    ASSERT_EQ(conn.cre(0), 5);
    ASSERT_EQ(conn.cre(1), 9);
    ASSERT_EQ(conn.cre(2), 137);

    AntisymFermionOnvConnection aconn(ket, bra);
    ASSERT_EQ(aconn.ncom(), 6);
    ASSERT_EQ(aconn.com(0), 1);
    ASSERT_EQ(aconn.com(1), 4);
    ASSERT_EQ(aconn.com(2), 6);
    ASSERT_EQ(aconn.com(3), 11);
    ASSERT_EQ(aconn.com(4), 19);
    ASSERT_EQ(aconn.com(5), 138);

    ASSERT_FALSE(aconn.phase());
}

TEST(Connection, Phase) {

    SparseArrayFileReader<float> file_reader(
            defs::assets_root + "/parity_test/parity_8.txt",
            16ul, false, false);
    defs::inds inds(16);
    float value;

    buffered::FermionOnv bra(4);
    buffered::FermionOnv ket(4);
    buffered::FermionOnv work_det(4);
    AntisymFermionOnvConnection connection(ket);

    while (file_reader.next(inds, value)) {
        bra.zero();
        ket.zero();
        for (size_t i = 0ul; i < 8ul; ++i) {
            if (inds[i]) bra.set(i);
        }
        for (size_t i = 8ul; i < 16ul; ++i) {
            if (inds[i]) ket.set(i - 8);
        }
        if (bra.is_zero() || ket.is_zero()) continue;
        if (bra.nsetbit() != ket.nsetbit()) continue;
        connection.connect(ket, bra);
        ASSERT_EQ(connection.phase(), value < 0);
        connection.apply(ket, work_det);
        ASSERT_TRUE(bra == work_det);
        ASSERT_EQ(connection.phase(), value < 0);
    }
}

TEST(Connection, MultiWordPhase) {
    // choose a large enough nsite so that multiple 64bit datawords are required
    const size_t nsite = 100;
    buffered::FermionOnv bra(nsite);
    buffered::FermionOnv ket(nsite);
    buffered::FermionOnv check_onv(nsite);
    AntisymFermionOnvConnection connection(nsite);

    bra = {{3,  4,  8,  13, 89},
           {13, 78, 95, 98, 99}};
    ket = bra;
    connection.connect(bra, ket);
    ASSERT_FALSE(connection.phase());
    connection.apply(bra, check_onv);
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //                              78 -> 3 (-1)         (-1)
    ket = {{3, 4,  8,  13, 89},
           {3, 13, 95, 98, 99}};
    connection.connect(bra, ket);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), nsite + 78);
    ASSERT_EQ(connection.cre(0), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());

    ASSERT_EQ(ket, check_onv);

    //                              98 -> 3 (-1)         (-1)
    ket = {{3, 4,  8,  13, 89},
           {3, 13, 78, 95, 99}};
    connection.connect(bra, ket);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), nsite + 98);
    ASSERT_EQ(connection.cre(0), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //     4 -> 99 (-1)             98 -> 3 (-1)         (+1)
    ket = {{3, 8,  13, 89, 99},
           {3, 13, 78, 95, 99}};
    connection.connect(bra, ket);
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(connection.ann(0), 4);
    ASSERT_EQ(connection.ann(1), nsite + 98);
    ASSERT_EQ(connection.cre(0), 99);
    ASSERT_EQ(connection.cre(1), nsite + 3);
    connection.apply(bra, check_onv);
    ASSERT_FALSE(connection.phase());
    ASSERT_EQ(ket, check_onv);

    //     4 -> 99 (-1)             98 -> 14 (+1)        (-1)
    ket = {{3,  8,  13, 89, 99},
           {13, 14, 78, 95, 99}};
    connection.connect(bra, ket);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(connection.ann(0), 4);
    ASSERT_EQ(connection.ann(1), nsite + 98);
    ASSERT_EQ(connection.cre(0), 99);
    ASSERT_EQ(connection.cre(1), nsite + 14);
    connection.apply(bra, check_onv);
    ASSERT_TRUE(connection.phase());
    ASSERT_EQ(ket, check_onv);
}
