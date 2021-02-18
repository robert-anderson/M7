//
// Created by rja on 06/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/nd/NdFormat.h"

TEST(NdArrayFormat, Test1D) {
    NdFormat<1> format(9);
    size_t iflat = 0ul;
    std::array<size_t, 1> iarr;
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        ASSERT_EQ(format.flatten(i0), iflat);
        format.decode_flat(iflat, iarr);
        ASSERT_EQ(iarr[0], i0);
        ++iflat;
    }
}

TEST(NdArrayFormat, Test2D) {
    NdFormat<2> format({4, 5});
    size_t iflat = 0ul;
    std::array<size_t, 2> iarr;
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            ASSERT_EQ(format.flatten(i0, i1), iflat);
            format.decode_flat(iflat, iarr);
            ASSERT_EQ(iarr[0], i0);
            ASSERT_EQ(iarr[1], i1);
            ++iflat;
        }
    }
}

TEST(NdArrayFormat, Test2DEqualExtents) {
    const size_t n = 7;
    NdFormat<2> format(n);
    size_t iflat = 0ul;
    std::array<size_t, 2> iarr;
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            ASSERT_EQ(format.flatten(i0, i1), iflat);
            format.decode_flat(iflat, iarr);
            ASSERT_EQ(iarr[0], i0);
            ASSERT_EQ(iarr[1], i1);
            ++iflat;
        }
    }
}

TEST(NdArrayFormat, Test3D) {
    NdFormat<3> format({4, 5, 3});
    size_t iflat = 0ul;
    std::array<size_t, 3> iarr;
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            for (size_t i2 = 0ul; i2 < format.extent(2); ++i2) {
                ASSERT_EQ(format.flatten(i0, i1, i2), iflat);
                format.decode_flat(iflat, iarr);
                ASSERT_EQ(iarr[0], i0);
                ASSERT_EQ(iarr[1], i1);
                ASSERT_EQ(iarr[2], i2);
                ++iflat;
            }
        }
    }
}

TEST(NdArrayFormat, Test3DEqualExtents) {
    const size_t n = 7;
    NdFormat<3> format(n);
    size_t iflat = 0ul;
    std::array<size_t, 3> iarr;
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            for (size_t i2 = 0ul; i2 < n; ++i2) {
                ASSERT_EQ(format.flatten(i0, i1, i2), iflat);
                format.decode_flat(iflat, iarr);
                ASSERT_EQ(iarr[0], i0);
                ASSERT_EQ(iarr[1], i1);
                ASSERT_EQ(iarr[2], i2);
                ++iflat;
            }
        }
    }
}

TEST(NdArrayFormat, Test4D) {
    NdFormat<4> format({4, 5, 3, 2});
    size_t iflat = 0ul;
    std::array<size_t, 4> iarr;
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            for (size_t i2 = 0ul; i2 < format.extent(2); ++i2) {
                for (size_t i3 = 0ul; i3 < format.extent(3); ++i3) {
                    ASSERT_EQ(format.flatten(i0, i1, i2, i3), iflat);
                    format.decode_flat(iflat, iarr);
                    ASSERT_EQ(iarr[0], i0);
                    ASSERT_EQ(iarr[1], i1);
                    ASSERT_EQ(iarr[2], i2);
                    ASSERT_EQ(iarr[3], i3);
                    ++iflat;
                }
            }
        }
    }
}

TEST(NdArrayFormat, Test4DEqualExtents) {
    const size_t n = 5;
    NdFormat<4> format(n);
    size_t iflat = 0ul;
    std::array<size_t, 4> iarr;
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            for (size_t i2 = 0ul; i2 < n; ++i2) {
                for (size_t i3 = 0ul; i3 < n; ++i3) {
                    ASSERT_EQ(format.flatten(i0, i1, i2, i3), iflat);
                    format.decode_flat(iflat, iarr);
                    ASSERT_EQ(iarr[0], i0);
                    ASSERT_EQ(iarr[1], i1);
                    ASSERT_EQ(iarr[2], i2);
                    ASSERT_EQ(iarr[3], i3);
                    ++iflat;
                }
            }
        }
    }
}
