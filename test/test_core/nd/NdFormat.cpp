//
// Created by Robert J. Anderson on 06/10/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/nd/NdFormat.h"

TEST(NdFormat, SubFormats) {
    NdFormat<4> format({3, 4, 2, 6});
    ASSERT_EQ(format.m_nelement, 3*4*2*6);
    ASSERT_EQ(format.major_dims().m_nelement, 3 * 4 * 2);
    ASSERT_EQ(format.major_dims().major_dims().m_nelement, 3 * 4);
    ASSERT_EQ(format.major_dims().major_dims().major_dims().m_nelement, 3);
    ASSERT_EQ(format.major_dims().major_dims().major_dims().major_dims().m_nelement, 1);
    ASSERT_EQ(format.minor_dims().m_nelement, 4 * 2 * 6);
    ASSERT_EQ(format.minor_dims().minor_dims().m_nelement, 2 * 6);
    ASSERT_EQ(format.minor_dims().minor_dims().minor_dims().m_nelement, 6);
    ASSERT_EQ(format.minor_dims().minor_dims().minor_dims().minor_dims().m_nelement, 1);

    ASSERT_EQ(format.major_dims<3>().m_nelement, 3 * 4 * 2);
    ASSERT_EQ(format.major_dims<2>().m_nelement, 3 * 4);
    ASSERT_EQ(format.major_dims<1>().m_nelement, 3);
    ASSERT_EQ(format.major_dims<0>().m_nelement, 1);
    ASSERT_EQ(format.minor_dims<3>().m_nelement, 4 * 2 * 6);
    ASSERT_EQ(format.minor_dims<2>().m_nelement, 2 * 6);
    ASSERT_EQ(format.minor_dims<1>().m_nelement, 6);
    ASSERT_EQ(format.minor_dims<0>().m_nelement, 1);

    auto major = format.major_dims<2>();
    auto minor = format.minor_dims<2>();

    uint_t iflat = 0ul;
    for (uint_t imajor_flat = 0ul; imajor_flat<major.m_nelement; ++imajor_flat){
        for (uint_t iminor_flat = 0ul; iminor_flat<minor.m_nelement; ++iminor_flat){
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
    uint_t iflat = 0ul;
    uinta_t<1> iarr;
    for (uint_t i0 = 0ul; i0 < format.m_shape[0]; ++i0) {
        ASSERT_EQ(format.flatten(i0), iflat);
        format.decode_flat(iflat, iarr);
        ASSERT_EQ(iarr[0], i0);
        ++iflat;
    }
}

TEST(NdFormat, Test2D) {
    NdFormat<2> format({4, 5});
    uint_t iflat = 0ul;
    uinta_t<2> iarr;
    for (uint_t i0 = 0ul; i0 < format.m_shape[0]; ++i0) {
        for (uint_t i1 = 0ul; i1 < format.m_shape[1]; ++i1) {
            ASSERT_EQ(format.flatten(i0, i1), iflat);
            format.decode_flat(iflat, iarr);
            ASSERT_EQ(iarr[0], i0);
            ASSERT_EQ(iarr[1], i1);
            ++iflat;
        }
    }
}

TEST(NdFormat, Test2DEqualExtents) {
    const uint_t n = 7;
    NdFormat<2> format(n);
    uint_t iflat = 0ul;
    uinta_t<2> iarr;
    for (uint_t i0 = 0ul; i0 < n; ++i0) {
        for (uint_t i1 = 0ul; i1 < n; ++i1) {
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
    uint_t iflat = 0ul;
    uinta_t<3> iarr;
    for (uint_t i0 = 0ul; i0 < format.m_shape[0]; ++i0) {
        for (uint_t i1 = 0ul; i1 < format.m_shape[1]; ++i1) {
            for (uint_t i2 = 0ul; i2 < format.m_shape[2]; ++i2) {
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
    const uint_t n = 7;
    NdFormat<3> format(n);
    uint_t iflat = 0ul;
    uinta_t<3> iarr;
    for (uint_t i0 = 0ul; i0 < n; ++i0) {
        for (uint_t i1 = 0ul; i1 < n; ++i1) {
            for (uint_t i2 = 0ul; i2 < n; ++i2) {
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
    uint_t iflat = 0ul;
    uinta_t<4> iarr;
    for (uint_t i0 = 0ul; i0 < format.m_shape[0]; ++i0) {
        for (uint_t i1 = 0ul; i1 < format.m_shape[1]; ++i1) {
            for (uint_t i2 = 0ul; i2 < format.m_shape[2]; ++i2) {
                for (uint_t i3 = 0ul; i3 < format.m_shape[3]; ++i3) {
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
    const uint_t n = 5;
    NdFormat<4> format(n);
    uint_t iflat = 0ul;
    uinta_t<4> iarr{};
    for (uint_t i0 = 0ul; i0 < n; ++i0) {
        for (uint_t i1 = 0ul; i1 < n; ++i1) {
            for (uint_t i2 = 0ul; i2 < n; ++i2) {
                for (uint_t i3 = 0ul; i3 < n; ++i3) {
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

TEST(NdFormat, CombineIndsFromSubformat31) {
    NdFormat<4> format({3, 4, 2, 6});
    auto major = format.major_dims<3>();
    auto minor = format.minor_dims<1>();
    ASSERT_EQ(major.m_nelement, 3*4*2);
    ASSERT_EQ(minor.m_nelement, 6);

    uint_t i = 0ul;
    for (uint_t imajor = 0ul; imajor<major.m_nelement; ++imajor){
        for (uint_t iminor = 0ul; iminor<minor.m_nelement; ++iminor) {
            ASSERT_EQ(format.combine<3>(imajor, iminor), i);
            ++i;
        }
    }
    ASSERT_EQ(i, format.m_nelement);
}

TEST(NdFormat, CombineIndsFromSubformat22) {
    NdFormat<4> format({3, 4, 2, 6});
    auto major = format.major_dims<2>();
    auto minor = format.minor_dims<2>();
    ASSERT_EQ(major.m_nelement, 3*4);
    ASSERT_EQ(minor.m_nelement, 2*6);

    uint_t i = 0ul;
    for (uint_t imajor = 0ul; imajor<major.m_nelement; ++imajor){
        for (uint_t iminor = 0ul; iminor<minor.m_nelement; ++iminor) {
            ASSERT_EQ(format.combine<2>(imajor, iminor), i);
            ++i;
        }
    }
    ASSERT_EQ(i, format.m_nelement);
}