//
// Created by rja on 06/06/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/FcidumpFileReader.h"

TEST(FcidumpFileReader, Real_6orb){
    FcidumpFileReader file_reader(defs::assets_root+"/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_FALSE(file_reader.m_spin_resolved);
    ASSERT_TRUE(file_reader.spin_conserving());
    ASSERT_EQ(file_reader.m_nspatorb, 6);
    defs::inds inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::inds test_inds(4);
    // first entry
    test_inds = {0,0,0,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v),  0.5406487462037872));
    // scan to arbitrary element
    for (size_t i=0; i<17; ++i) file_reader.next(inds, v);
    test_inds = {2,1,4,3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), 0.01759459248922075));
    // scan to final element
    while(file_reader.next(inds, v)){}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), -98.46644393370157));
}

TEST(FcidumpFileReader, Integer_8orb){
    FcidumpFileReader file_reader(defs::assets_root+"/Hubbard_U4_8site/FCIDUMP", false);
    ASSERT_FALSE(file_reader.m_spin_resolved);
    ASSERT_TRUE(file_reader.spin_conserving());
    ASSERT_EQ(file_reader.m_nspatorb, 8);
    defs::inds inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::inds test_inds(4);
    // core energy is the first entry
    test_inds = {~0ul,~0ul,~0ul,~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), 0.0));
    // scan to arbitrary element
    for (size_t i=0; i<8; ++i) file_reader.next(inds, v);
    test_inds = {7,7,7,7};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), 4.0));
    // scan to final element
    while(file_reader.next(inds, v)){}
    test_inds = {7, 6, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), -1.0));
}

#ifdef ENABLE_COMPLEX
TEST(FcidumpFileReader, Complex_10orb){
    FcidumpFileReader file_reader(defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_EQ(file_reader.m_nspatorb, 5);
    defs::inds inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::inds test_inds(4);
    // first entry
    test_inds = {0,0,0,0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), 2.2752637995109302));
    ASSERT_TRUE(consts::floats_equal(consts::imag(v), 0.0));
    // scan to arbitrary element
    for (size_t i=0; i<20; ++i) file_reader.next(inds, v);
    // (-0.00851916802083687,-0.005287130898791)   5   3   7   1
    test_inds = {5,3,7,1};
    file_reader.inds_to_orbs(test_inds);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(v), -0.00851916802083687));
    ASSERT_TRUE(consts::floats_equal(consts::imag(v), -0.005287130898791));
    // scan to final element
    while(file_reader.next(inds, v)){}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::float_is_zero(v));
}
#endif