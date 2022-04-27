//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "M7_lib/io/EbdumpFileReader.h"
#include "gtest/gtest.h"

TEST(EbDumpFileReader, ReadFile){
    std::string fname = defs::assets_root + "/Hubbard_U4_4site/EBDUMP";
    EbdumpFileReader file_reader(fname);
    ASSERT_FALSE(file_reader.m_header.m_uhf);
    ASSERT_EQ(file_reader.m_header.m_nmode, 4ul);
    ASSERT_EQ(file_reader.m_header.m_nsite, 4ul);
    defs::inds inds(3);
    defs::inds test_inds(3);
    defs::ham_t value;

    test_inds = {0, 0, 0};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), -0.01832284815866287);

    // scan to arbitrary element
    for (size_t i=0; i<18; ++i) file_reader.next(inds, value);
    test_inds = {1, 0, 2};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), -0.009244953059296356);

    // scan to final element
    while(file_reader.next(inds, value)){}
    test_inds = {3, 3, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), -0.09669569380966529);
}

TEST(EbDumpFileReader, SpinResolved){
    std::string fname = defs::assets_root + "/SpinResolvedEbdump/EBDUMP";
    EbdumpFileReader file_reader(fname);
    ASSERT_TRUE(file_reader.m_header.m_uhf);
    ASSERT_EQ(file_reader.m_header.m_nmode, 6ul);
    ASSERT_EQ(file_reader.m_header.m_nsite, 6ul);

    defs::inds inds(3);
    defs::inds test_inds(3);
    defs::ham_t value;

    test_inds = {0, 0, 6};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), 0.07953916361235415);

#if 0
    // scan to arbitrary element
    for (size_t i=0; i<18; ++i) file_reader.next(inds, value);
    test_inds = {1, 0, 2};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), -0.009244953059296356);

    // scan to final element
    while(file_reader.next(inds, value)){}
    test_inds = {3, 3, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(consts::real(value), -0.09669569380966529);
#endif
}