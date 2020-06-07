//
// Created by rja on 06/06/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/FcidumpFileIterator.h"

TEST(FcidumpFileIterator, Real_6orb){
    typedef double T;
    FcidumpFileIterator<T> file_iterator(defs::assets_root+"/RHF_N2_6o6e/FCIDUMP");
    ASSERT_EQ(file_iterator.nsite(), 6);
    defs::inds inds(4);
    T v;
    file_iterator.next(inds, v);
    defs::inds test_inds(4);
    // first entry
    test_inds = {0,0,0,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(v, 0.53644984080846303));
    // scan to arbitrary element
    for (size_t i=0; i<17; ++i) file_iterator.next(inds, v);
    test_inds = {4,5,2,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(v, 0.0043863275703943001));
    // scan to final element
    while(file_iterator.next(inds, v)){}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(v, -98.3339467443989));
}


TEST(FcidumpFileIterator, Complex_10orb){
    typedef std::complex<double> T;
    FcidumpFileIterator<T> file_iterator(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP");
    ASSERT_EQ(file_iterator.nsite(), 5);
    defs::inds inds(4);
    T v;
    file_iterator.next(inds, v);
    defs::inds test_inds(4);
    // first entry
    test_inds = {0,0,0,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(v, T(2.2752637995109302, 2.379e-17)));
    // scan to arbitrary element
    for (size_t i=0; i<20; ++i) file_iterator.next(inds, v);
    test_inds = {4,2,6,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(v, T(-0.00851916802083687,-0.005287130898791)));
    // scan to final element
    while(file_iterator.next(inds, v)){}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::float_is_zero(v));
}