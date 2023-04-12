//
// Created by rja on 04/06/22.
//

#include "test_core/defs.h"
#include "M7_lib/integrals/IntegralArray1e.h"


TEST(IntegralArray1e, SymNone_real) {
    typedef double T;
    integrals_1e::SymNone<T> ints(6);
    if (mpi::on_node_i_am_root()) {
        ASSERT_TRUE(ints.set_(2, 5, 0.123));
        ASSERT_TRUE(ints.set_(2, 5, 0.123));
        ASSERT_FALSE(ints.set_(2, 5, 0.1234));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(2, 5), 0.123);
    ASSERT_NEAR_EQ(ints.get(5, 2), 0.0);
    mpi::barrier_on_node();
    if (mpi::on_node_i_am_root()) {
        ASSERT_TRUE(ints.set_(5, 2, 0.234));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(5, 2), 0.234);
}

TEST(IntegralArray1e, SymNone_complex) {
    typedef std::complex<double> T;
    integrals_1e::SymNone<T> ints(6);
    const T v = {0.123, -0.234};
    if (mpi::on_node_i_am_root()) {
        ASSERT_TRUE(ints.set_(2, 5, v));
        ASSERT_TRUE(ints.set_(2, 5, v));
        ASSERT_FALSE(ints.set_(2, 5, {0.1234, 0.2345}));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(2, 5), v);
    ASSERT_NEAR_EQ(ints.get(5, 2), 0.0);
}

TEST(IntegralArray1e, SymH_real) {
    typedef double T;
    integrals_1e::SymH<T> ints(6);
    if (mpi::on_node_i_am_root()) {
        ASSERT_TRUE(ints.set_(2, 5, 0.123));
        ASSERT_TRUE(ints.set_(2, 5, 0.123));
        ASSERT_FALSE(ints.set_(2, 5, 0.1234));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(2, 5), 0.123);
    ASSERT_NEAR_EQ(ints.get(5, 2), 0.123);
    mpi::barrier_on_node();
    if (mpi::on_node_i_am_root()) {
        ASSERT_FALSE(ints.set_(5, 2, 0.234));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(5, 2), 0.123);
}

TEST(IntegralArray1e, SymH_complex) {
    typedef std::complex<double> T;
    integrals_1e::SymH<T> ints(6);
    const T v = {0.123, -0.234};
    const T vc = {0.123, 0.234};
    if (mpi::on_node_i_am_root()) {
        ASSERT_TRUE(ints.set_(2, 5, v));
        ASSERT_TRUE(ints.set_(2, 5, v));
        ASSERT_FALSE(ints.set_(2, 5, {0.1234, 0.2345}));
    }
    mpi::barrier_on_node();
    ASSERT_NEAR_EQ(ints.get(2, 5), v);
    ASSERT_NEAR_EQ(ints.get(5, 2), vc);
    mpi::barrier_on_node();
    ASSERT_FALSE(ints.set_(5, 2, v));
    ASSERT_TRUE(ints.set_(5, 2, vc));
}