//
// Created by Robert J. Anderson on 06/06/2020.
//

#include "test_core/defs.h"
#include "M7_lib/io/FcidumpFileReader.h"

TEST(FcidumpInfo, EmptyFilename) {
    FcidumpInfo header("");
    ASSERT_EQ(header.m_nsite, 0ul);
    ASSERT_EQ(header.m_nelec, 0ul);
    ASSERT_EQ(header.m_spin_resolved, false);
    ASSERT_EQ(header.m_uhf, false);
    ASSERT_EQ(header.m_orbsym.size(), 0ul);
    ASSERT_EQ(header.m_relativistic, false);
}

TEST(FcidumpFileReader, Real_6orb) {
    FcidumpFileReader file_reader(PROJECT_ROOT"/assets/RHF_N2_6o6e/FCIDUMP", false);
    ASSERT_FALSE(file_reader.m_info.m_spin_resolved);
    ASSERT_TRUE(file_reader.spin_conserving());
    ASSERT_EQ(file_reader.m_info.m_nsite, 6);
    defs::ivec_t inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::ivec_t test_inds(4);
    // first entry
    test_inds = {0, 0, 0, 0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 0.5406487462037872);
    // scan to arbitrary element
    for (size_t i = 0; i < 17; ++i) file_reader.next(inds, v);
    test_inds = {2, 1, 4, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 0.01759459248922075);
    // scan to final element
    while (file_reader.next(inds, v)) {}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, -98.46644393370157);
}

TEST(FcidumpFileReader, Integer_8orb) {
    FcidumpFileReader file_reader(PROJECT_ROOT"/assets/Hubbard_U4_8site/FCIDUMP", false);
    ASSERT_FALSE(file_reader.m_info.m_spin_resolved);
    ASSERT_TRUE(file_reader.spin_conserving());
    ASSERT_EQ(file_reader.m_info.m_nsite, 8);
    defs::ivec_t inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::ivec_t test_inds(4);
    // core energy is the first entry
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 0.0);
    // scan to arbitrary element
    for (size_t i = 0; i < 8; ++i) file_reader.next(inds, v);
    test_inds = {7, 7, 7, 7};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 4.0);
    // scan to final element
    while (file_reader.next(inds, v)) {}
    test_inds = {7, 6, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, -1.0);
}

TEST(FcidumpFileReader, Molcas) {
    FcidumpFileReader file_reader(PROJECT_ROOT"/assets/O2_Molcas/FCIDUMP", false);
    ASSERT_FALSE(file_reader.m_info.m_spin_resolved);
    ASSERT_TRUE(file_reader.spin_conserving());
    ASSERT_EQ(file_reader.m_info.m_nsite, 6);
    defs::ivec_t inds(4);
    defs::ham_t v;
    file_reader.next(inds, v);
    defs::ivec_t test_inds(4);

    test_inds = {0, 0, 0, 0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 0.75132124044);
    file_reader.next(inds, v);
    test_inds = {1, 0, 1, 0};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(v, 0.32107592937E-01);
}

#ifdef ENABLE_COMPLEX
TEST(FcidumpFileReader, Complex_10orb){
    FcidumpFileReader file_reader(PROJECT_ROOT"/assets/DHF_Be_STO-3G/FCIDUMP", false);
    ASSERT_EQ(file_reader.m_nspatorb, 5);
    defs::ivec_t ivec_t(4);
    defs::ham_t v;
    file_reader.next(ivec_t, v);
    defs::ivec_t test_inds(4);
    // first entry
    test_inds = {0,0,0,0};
    ASSERT_TRUE(std::equal(ivec_t.begin(), ivec_t.end(), test_inds.begin()));
    ASSERT_TRUE(datatype::floats_equal(datatype::real(v), 2.2752637995109302));
    ASSERT_TRUE(datatype::floats_equal(datatype::imag(v), 0.0));
    // scan to arbitrary element
    for (size_t i=0; i<20; ++i) file_reader.next(ivec_t, v);
    // (-0.00851916802083687,-0.005287130898791)   5   3   7   1
    test_inds = {5,3,7,1};
    file_reader.inds_to_orbs(test_inds);
    ASSERT_TRUE(std::equal(ivec_t.begin(), ivec_t.end(), test_inds.begin()));
    ASSERT_TRUE(datatype::floats_equal(datatype::real(v), -0.00851916802083687));
    ASSERT_TRUE(datatype::floats_equal(datatype::imag(v), -0.005287130898791));
    // scan to final element
    while(file_reader.next(ivec_t, v)){}
    test_inds = {~0ul, ~0ul, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(ivec_t.begin(), ivec_t.end(), test_inds.begin()));
    ASSERT_TRUE(datatype::float_is_zero(v));
}
#endif
