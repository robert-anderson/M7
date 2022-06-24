//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "M7_lib/io/EbdumpFileReader.h"
#include "gtest/gtest.h"
#include "test_core/defs.h"

TEST(EbDumpFileReader, ReadFile){
    std::string fname = PROJECT_ROOT"/assets/Hubbard_U4_4site/EBDUMP";
    EbdumpFileReader file_reader(fname);
    ASSERT_FALSE(file_reader.m_info.m_uhf);
    ASSERT_EQ(file_reader.m_info.m_nmode, 4ul);
    ASSERT_EQ(file_reader.m_info.m_nsite, 4ul);
    uintv_t inds(3);
    uintv_t test_inds(3);
    ham_t value;

    test_inds = {0, 0, 0};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, -0.01832284815866287);

    // scan to arbitrary element
    for (uint_t i=0; i<18; ++i) file_reader.next(inds, value);
    test_inds = {1, 0, 2};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, -0.009244953059296356);

    // scan to final element
    while(file_reader.next(inds, value)){}
    test_inds = {3, 3, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, -0.09669569380966529);
}

TEST(EbDumpFileReader, SpinResolved){
    std::string fname = PROJECT_ROOT"/assets/SpinResolvedEbdump/EBDUMP";
    EbdumpFileReader file_reader(fname);
    ASSERT_TRUE(file_reader.m_info.m_uhf);
    ASSERT_EQ(file_reader.m_info.m_nmode, 6ul);
    ASSERT_EQ(file_reader.m_info.m_nsite, 6ul);

    uintv_t inds(3);
    uintv_t test_inds(3);
    ham_t value;

    // 1 1 7 -> 0 0 6 -> 0 0 3
    test_inds = {0, 0, 3};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, 0.07953916361235415);

    // scan to arbitrary element
    for (uint_t i=0; i<18; ++i) file_reader.next(inds, value);
    // 3 1 7 -> 2 0 6 -> 2 0 3
    test_inds = {2, 0, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, 0.1575692477527167);

    // scan to final element
    while(file_reader.next(inds, value)){}
    // 6 12 6 -> 5 11 5 -> 5 11 8
    test_inds = {5, 11, 8};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_NEARLY_EQ(value, -0.03197691666078238);
}