//
// Created by rja on 06/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/nd/NdArrayFormat.h"

TEST(NdArrayFormat, Test1D) {
    NdArrayFormat<1> format(9);
    size_t iflat = 0ul;
    defs::inds inds(1, 0ul);
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        ASSERT_EQ(format.flat(i0), iflat);
        format.decode_flat(iflat, inds);
        ASSERT_EQ(inds[0], i0);
        ++iflat;
    }
}

TEST(NdArrayFormat, Test2D) {
    NdArrayFormat<2> format(4, 5);
    size_t iflat = 0ul;
    defs::inds inds(2, 0ul);
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            ASSERT_EQ(format.flat(i0, i1), iflat);
            format.decode_flat(iflat, inds);
            ASSERT_EQ(inds[0], i0);
            ASSERT_EQ(inds[1], i1);
            ++iflat;
        }
    }
}

TEST(NdArrayFormat, Test2DEqualExtents) {
    const size_t n = 7;
    NdArrayFormat<2> format(n);
    size_t iflat = 0ul;
    defs::inds inds(2, 0ul);
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            ASSERT_EQ(format.flat(i0, i1), iflat);
            format.decode_flat(iflat, inds);
            ASSERT_EQ(inds[0], i0);
            ASSERT_EQ(inds[1], i1);
            ++iflat;
        }
    }
}

TEST(NdArrayFormat, Test3D) {
    NdArrayFormat<3> format(4, 5, 3);
    size_t iflat = 0ul;
    defs::inds inds(3, 0ul);
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            for (size_t i2 = 0ul; i2 < format.extent(2); ++i2) {
                ASSERT_EQ(format.flat(i0, i1, i2), iflat);
                format.decode_flat(iflat, inds);
                ASSERT_EQ(inds[0], i0);
                ASSERT_EQ(inds[1], i1);
                ASSERT_EQ(inds[2], i2);
                ++iflat;
            }
        }
    }
}

TEST(NdArrayFormat, Test3DEqualExtents) {
    const size_t n = 7;
    NdArrayFormat<3> format(n);
    size_t iflat = 0ul;
    defs::inds inds(3, 0ul);
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            for (size_t i2 = 0ul; i2 < n; ++i2) {
                ASSERT_EQ(format.flat(i0, i1, i2), iflat);
                format.decode_flat(iflat, inds);
                ASSERT_EQ(inds[0], i0);
                ASSERT_EQ(inds[1], i1);
                ASSERT_EQ(inds[2], i2);
                ++iflat;
            }
        }
    }
}

TEST(NdArrayFormat, Test4D) {
    NdArrayFormat<4> format(4, 5, 3, 2);
    size_t iflat = 0ul;
    defs::inds inds(4, 0ul);
    for (size_t i0 = 0ul; i0 < format.extent(0); ++i0) {
        for (size_t i1 = 0ul; i1 < format.extent(1); ++i1) {
            for (size_t i2 = 0ul; i2 < format.extent(2); ++i2) {
                for (size_t i3 = 0ul; i3 < format.extent(3); ++i3) {
                    ASSERT_EQ(format.flat(i0, i1, i2, i3), iflat);
                    format.decode_flat(iflat, inds);
                    ASSERT_EQ(inds[0], i0);
                    ASSERT_EQ(inds[1], i1);
                    ASSERT_EQ(inds[2], i2);
                    ASSERT_EQ(inds[3], i3);
                    ++iflat;
                }
            }
        }
    }
}

TEST(NdArrayFormat, Test4DEqualExtents) {
    const size_t n = 5;
    NdArrayFormat<4> format(n);
    size_t iflat = 0ul;
    defs::inds inds(4, 0ul);
    for (size_t i0 = 0ul; i0 < n; ++i0) {
        for (size_t i1 = 0ul; i1 < n; ++i1) {
            for (size_t i2 = 0ul; i2 < n; ++i2) {
                for (size_t i3 = 0ul; i3 < n; ++i3) {
                    ASSERT_EQ(format.flat(i0, i1, i2, i3), iflat);
                    format.decode_flat(iflat, inds);
                    ASSERT_EQ(inds[0], i0);
                    ASSERT_EQ(inds[1], i1);
                    ASSERT_EQ(inds[2], i2);
                    ASSERT_EQ(inds[3], i3);
                    ++iflat;
                }
            }
        }
    }
}
