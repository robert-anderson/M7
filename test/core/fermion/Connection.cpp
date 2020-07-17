//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/fermion/Determinant.h>
#include <src/core/fermion/Connection.h>
#include <src/core/io/SparseArrayFileReader.h>
#include "gtest/gtest.h"

TEST(Connection, ParticleNumberConserving){
    const size_t nsite = 20;
    Determinant ket(nsite);
    Determinant bra(nsite);

    defs::inds ketoccorbs = {1, 4, 6, 8, 11, 19, 20, 38, 39};

    defs::inds braoccorbs = {1, 4, 5, 6, 9, 11, 19, 37, 38};
    ket.set(ketoccorbs);
    bra.set(braoccorbs);

    ASSERT_EQ(ket.nsetbit(), 9);
    ASSERT_EQ(bra.nsetbit(), 9);

    Connection conn(ket, bra);

    ASSERT_EQ(conn.nann(), 3);
    ASSERT_EQ(conn.ncre(), 3);

    ASSERT_EQ(conn.ann(0), 8);
    ASSERT_EQ(conn.ann(1), 20);
    ASSERT_EQ(conn.ann(2), 39);

    ASSERT_EQ(conn.cre(0), 5);
    ASSERT_EQ(conn.cre(1), 9);
    ASSERT_EQ(conn.cre(2), 37);

    AntisymConnection aconn(ket, bra);
    ASSERT_EQ(aconn.ncom(), 6);
    ASSERT_EQ(aconn.com(0), 1);
    ASSERT_EQ(aconn.com(1), 4);
    ASSERT_EQ(aconn.com(2), 6);
    ASSERT_EQ(aconn.com(3), 11);
    ASSERT_EQ(aconn.com(4), 19);
    ASSERT_EQ(aconn.com(5), 38);

    ASSERT_FALSE(aconn.phase());
}

TEST(Connection, Phase) {

    SparseArrayFileReader<float> file_reader(defs::assets_root + "/parity_test/parity_8.txt", 16);
    defs::inds inds(16);
    float value;

    Determinant bra(4);
    Determinant ket(4);
    Determinant work_det(4);
    AntisymConnection connection(ket);

    while (file_reader.next(inds, value)) {
        bra.zero();
        ket.zero();
        for (auto i{0ul}; i < 8; ++i) {
            if (!inds[i]) bra.set(i);
        }
        for (auto i{8ul}; i < 16; ++i) {
            if (!inds[i]) ket.set(i - 8);
        }
        if (bra.is_zero() || ket.is_zero()) continue;
        if (bra.nsetbit() != ket.nsetbit()) continue;
        connection.connect(ket, bra);
        ASSERT_EQ(connection.phase(), value < 0);
        connection.apply(ket, work_det);
        ASSERT_EQ(connection.phase(), value < 0);
        ASSERT_TRUE(bra==work_det);
    }
}