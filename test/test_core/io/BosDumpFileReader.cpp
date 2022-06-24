//
// Created by Robert J. Anderson on 22/08/2021.
//

#include "M7_lib/io/BosdumpFileReader.h"
#include "gtest/gtest.h"

TEST(BosdumpFileReader, ReadFile){
    std::string fname = PROJECT_ROOT"/assets/LandauLevels_5_5_15/BOSDUMP";
    BosdumpFileReader file_reader(fname);
    ASSERT_EQ(file_reader.m_header.m_nmode, 5ul);
    uintv_t inds(4);
    uintv_t test_inds(4);
    ham_t value;

    test_inds = {0, 0, 0, 0};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(arith::real(value), 1.0);

    // scan to arbitrary element
    // 0.3535533906 3 2 1 2
    for (uint_t i=3; i<37; ++i) file_reader.next(inds, value);
    test_inds = {2, 1, 0, 1};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(arith::real(value), 0.3535533906);

    // scan to final element
    // 0.2734375000 5 5 5 5
    while(file_reader.next(inds, value)){}
    test_inds = {4, 4, 4, 4};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_FLOAT_EQ(arith::real(value), 0.2734375000);
}