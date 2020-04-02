//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/fermion/Determinant.h>
#include <src/core/fermion/Connection.h>
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

    ASSERT_EQ(conn.m_ndes, 3);
    ASSERT_EQ(conn.m_ncre, 3);

    ASSERT_EQ(conn.m_des[0], 8);
    ASSERT_EQ(conn.m_des[1], 20);
    ASSERT_EQ(conn.m_des[2], 39);

    ASSERT_EQ(conn.m_cre[0], 5);
    ASSERT_EQ(conn.m_cre[1], 9);
    ASSERT_EQ(conn.m_cre[2], 37);

    AntisymConnection aconn(ket, bra);
    ASSERT_EQ(aconn.m_ncom, 6);
    ASSERT_EQ(aconn.m_com[0], 1);
    ASSERT_EQ(aconn.m_com[1], 4);
    ASSERT_EQ(aconn.m_com[2], 6);
    ASSERT_EQ(aconn.m_com[3], 11);
    ASSERT_EQ(aconn.m_com[4], 19);
    ASSERT_EQ(aconn.m_com[5], 38);

    ASSERT_FALSE(aconn.m_phase);
}