//
// Created by rja on 22/08/2021.
//

#include "src/core/io/BosdumpFileReader.h"
#include "gtest/gtest.h"

TEST(BosdumpFileReader, ReadFile){
    std::string fname = defs::assets_root + "/Hubbard_U4_4site/BOSDUMP";
    BosdumpFileReader file_reader(fname);
    ASSERT_FALSE(file_reader.m_spin_resolved);
    ASSERT_EQ(file_reader.m_norb, 4ul);
    defs::inds inds(2);
    defs::inds test_inds(2);
    defs::ham_t value;

    test_inds = {0, 0};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), 0.7134575551055222));

    // scan to arbitrary element
    for (size_t i=0; i<2; ++i) file_reader.next(inds, value);
    test_inds = {2, 2};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), 0.7835923888931371));

    // scan to final element
    while(file_reader.next(inds, value)){}
    test_inds = {3, 3};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), 0.708348137895966));
}