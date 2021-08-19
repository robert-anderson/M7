//
// Created by rja on 19/08/2021.
//

#include "src/core/io/EbdumpFileReader.h"
#include "gtest/gtest.h"

TEST(EbdumpFileReader, ReadFile){
    std::string fname = "/home/rja/CLionProjects/M7/assets/Hubbard_U4_4site/EBDUMP";
    EbdumpFileReader file_reader(fname);
    ASSERT_FALSE(file_reader.m_spin_resolved);
    ASSERT_EQ(file_reader.m_norb, 4ul);
    defs::inds inds(3);
    defs::inds test_inds(3);
    defs::ham_t value;

    test_inds = {0, 0, 0};
    file_reader.next(inds, value);
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), -0.01832284815866287));

    // scan to arbitrary element
    for (size_t i=0; i<18; ++i) file_reader.next(inds, value);
    test_inds = {1, 0, 2};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), -0.009244953059296356));

    // scan to first "omega" element
    while(file_reader.next(inds, value)){if (HamiltonianFileReader::nset_ind(inds)==1) break;}

    test_inds = {0, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), 0.7134575551055222));

    // scan to final element
    while(file_reader.next(inds, value)){}
    test_inds = {3, ~0ul, ~0ul};
    ASSERT_TRUE(std::equal(inds.begin(), inds.end(), test_inds.begin()));
    ASSERT_TRUE(consts::floats_equal(consts::real(value), 0.708348137895966));
}