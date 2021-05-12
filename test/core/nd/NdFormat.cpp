//
// Created by rja on 06/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/nd/NdFormat.h"

TEST(NdFormat, SubFormats) {
    NdFormat<4> format({3, 4, 2, 6});
    ASSERT_EQ(format.nelement(), 3*4*2*6);
    ASSERT_EQ(format.major_dims().nelement(), 3 * 4 * 2);
    ASSERT_EQ(format.major_dims().major_dims().nelement(), 3 * 4);
    ASSERT_EQ(format.major_dims().major_dims().major_dims().nelement(), 3);
    ASSERT_EQ(format.major_dims().major_dims().major_dims().major_dims().nelement(), 1);
    ASSERT_EQ(format.minor_dims().nelement(), 4 * 2 * 6);
    ASSERT_EQ(format.minor_dims().minor_dims().nelement(), 2 * 6);
    ASSERT_EQ(format.minor_dims().minor_dims().minor_dims().nelement(), 6);
    ASSERT_EQ(format.minor_dims().minor_dims().minor_dims().minor_dims().nelement(), 1);

    ASSERT_EQ(format.major_dims<1>().nelement(), 3 * 4 * 2);
    ASSERT_EQ(format.major_dims<2>().nelement(), 3 * 4);
    ASSERT_EQ(format.major_dims<3>().nelement(), 3);
    ASSERT_EQ(format.major_dims<4>().nelement(), 1);
    ASSERT_EQ(format.minor_dims<1>().nelement(), 4 * 2 * 6);
    ASSERT_EQ(format.minor_dims<2>().nelement(), 2 * 6);
    ASSERT_EQ(format.minor_dims<3>().nelement(), 6);
    ASSERT_EQ(format.minor_dims<4>().nelement(), 1);

    auto major = format.major_dims<2>();
    auto minor = format.minor_dims<2>();

    size_t iflat = 0ul;
    for (size_t imajor_flat = 0ul; imajor_flat<major.nelement(); ++imajor_flat){
        for (size_t iminor_flat = 0ul; iminor_flat<minor.nelement(); ++iminor_flat){
            ASSERT_EQ(iflat++, format.flatten<2>(imajor_flat, iminor_flat));
        }
    }
}

TEST(NdFormat, NamedDimensions) {
    NdFormat<3> format({3, 4, 2},
                       {"first", "second", "third"});
    ASSERT_EQ(format.to_string(), "first (3) second (4) third (2) ");
}

TEST(NdFormat, Test1D) {
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

TEST(NdFormat, Test2D) {
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

TEST(NdFormat, Test2DEqualExtents) {
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

TEST(NdFormat, Test3D) {
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

TEST(NdFormat, Test3DEqualExtents) {
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

TEST(NdFormat, Test4D) {
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

TEST(NdFormat, Test4DEqualExtents) {
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
