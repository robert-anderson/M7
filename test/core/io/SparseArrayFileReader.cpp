//
// Created by rja on 10/07/2020.
//

#include <src/core/util/defs.h>
#include "gtest/gtest.h"
#include <fstream>
#include <algorithm>

#include "src/core/io/FcidumpFileReader.h"

TEST(SparseArrayFileReader, IndexReading) {
    std::string line = "1 3 5 7";
    const char *p = line.begin().base();
    for (int i = 0; i < 4; ++i) {
        ASSERT_EQ(FcidumpFileReader<float>::read_unsigned(p), 2 * i + 1);
    }
    ASSERT_EQ(FcidumpFileReader<float>::read_unsigned(p), ~0ul);
}

TEST(SparseArrayFileReader, EntryReadingComplexIndsLast) {
    SparseArrayFileReader<std::complex<double>> file_reader(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP", 4);
    ASSERT_TRUE(file_reader.m_complex_valued);
    defs::inds inds(4);
    std::complex<double> v;
    file_reader.next(inds, v);
    ASSERT_EQ(inds, defs::inds({1, 1, 1, 1}));
    ASSERT_EQ(v, std::complex<double>(2.27526379951093, 2.37934998859424e-17));
    // on source file line #5 - advance 600 lines to #605
    // (-0.00231725227868085,-0.00527081024253388)   5   4   1   8
    for (size_t i = 0; i < 600; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({5, 4, 1, 8}));
    ASSERT_EQ(v, std::complex<double>(-0.00231725227868085, -0.00527081024253388));
    // last line is #2306
    for (size_t i = 0; i < 2306 - 605; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({0, 0, 0, 0}));
    ASSERT_EQ(v, std::complex<double>(0.0, 0.0));
    ASSERT_FALSE(file_reader.next(inds, v));
}

TEST(SparseArrayFileReader, EntryReadingRealIndsLastToComplexContainer) {
    FcidumpFileReader<std::complex<double>> file_reader(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", 4);
    defs::inds inds(4);
    std::complex<double> v;
    file_reader.next(inds, v);
    ASSERT_EQ(inds, defs::inds({1, 1, 1, 1}));
    ASSERT_EQ(v, std::complex<double>(0.536449840808463, 0));
    // on source file line #5 - advance 250 lines to #255
    //    0.0287053880352757   6   4   4   6
    for (size_t i = 0; i < 250; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({6, 4, 4, 6}));
    ASSERT_EQ(v, std::complex<double>(0.0287053880352757, 0.0));
    // last line is #275
    for (size_t i = 0; i < 275 - 255; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({0, 0, 0, 0}));
    ASSERT_EQ(v, std::complex<double>(-98.3339467443989, 0.0));
    ASSERT_FALSE(file_reader.next(inds, v));
}

TEST(SparseArrayFileReader, EntryReadingRealIndsLastToRealContainer) {
    SparseArrayFileReader<double> file_reader(defs::assets_root + "/RHF_N2_6o6e/FCIDUMP", 4);
    defs::inds inds(4);
    double v;
    file_reader.next(inds, v);
    ASSERT_EQ(inds, defs::inds({1, 1, 1, 1}));
    ASSERT_EQ(v, 0.536449840808463);
    // on source file line #5 - advance 250 lines to #255
    //    0.0287053880352757   6   4   4   6
    for (size_t i = 0; i < 250; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({6, 4, 4, 6}));
    ASSERT_EQ(v, 0.0287053880352757);
    // last line is #275
    for (size_t i = 0; i < 275 - 255; i++)
        ASSERT_TRUE(file_reader.next(inds, v));
    ASSERT_EQ(inds, defs::inds({0, 0, 0, 0}));
    ASSERT_EQ(v, -98.3339467443989);
    ASSERT_FALSE(file_reader.next(inds, v));
}