//
// Created by rja on 09/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

TEST(BufferedFields, CopyingNumbers){
    std::array<size_t, 2> shape {2, 3};
    buffered::Numbers<size_t, 2> field(shape);
    for (size_t i=0ul; i<field.m_nelement; ++i) field[i] = hashing::in_range(i, 10, 20);
    auto field_cpy = field;
    ASSERT_EQ(field_cpy.m_format.m_shape, shape);
    for (size_t i=0ul; i<field_cpy.m_nelement; ++i)
        ASSERT_EQ(field_cpy[i], hashing::in_range(i, 10, 20));
    ASSERT_EQ(field, field_cpy);
}

TEST(BufferedFields, MovingNumbers){
    std::array<size_t, 2> shape {2, 3};
    buffered::Numbers<size_t, 2> field(shape);
    for (size_t i=0ul; i<field.m_nelement; ++i) field[i] = hashing::in_range(i, 10, 20);
    auto field_mvd = std::move(field);
    ASSERT_EQ(field_mvd.m_format.m_shape, shape);
    for (size_t i=0ul; i<field_mvd.m_nelement; ++i)
        ASSERT_EQ(field_mvd[i], hashing::in_range(i, 10, 20));
}

TEST(BufferedFields, CopyingNdBitset){
    std::array<size_t, 3> shape {4, 3, 2};
    buffered::NdBitset<size_t, 3> field(shape);
    const size_t nsetbit = 8;
    ASSERT_EQ(field.m_format.m_nelement, 4*3*2);
    for (size_t i=0ul; i<nsetbit; ++i) field.set(hashing::in_range(i, 0, field.m_format.m_nelement));
    auto field_cpy = field;
    ASSERT_EQ(field_cpy.m_format.m_shape, shape);
    for (size_t i=0ul; i<nsetbit; ++i)
        ASSERT_TRUE(field_cpy.get(hashing::in_range(i, 0, field_cpy.m_format.m_nelement)));
    ASSERT_EQ(field, field_cpy);
}

TEST(BufferedFields, MovingNdBitset){
    std::array<size_t, 3> shape {4, 3, 2};
    buffered::NdBitset<size_t, 3> field(shape);
    const size_t nsetbit = 8;
    ASSERT_EQ(field.m_format.m_nelement, 4*3*2);
    for (size_t i=0ul; i<nsetbit; ++i) field.set(hashing::in_range(i, 0, field.m_format.m_nelement));
    auto field_mvd = std::move(field);
    ASSERT_EQ(field_mvd.m_format.m_shape, shape);
    for (size_t i=0ul; i<nsetbit; ++i)
        ASSERT_TRUE(field_mvd.get(hashing::in_range(i, 0, field_mvd.m_format.m_nelement)));
}