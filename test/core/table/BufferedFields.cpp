//
// Created by rja on 09/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"


/*
 * Implemented buffered fields of T (B<T>, where T is any class derived from FieldBase or CompositeField) must be:
 *  1. constructable via argument forwarding to T::T
 *  2. compatible with a call to any public method of T without cast to T&
 *  3. copy-constructable from const B<T>&
 *  4. copy-constructable from const T&
 *  5. copy-assignable from const B<T>&
 *  6. copy-assignable from const T&
 *  7. assignable via any other method T::operator= without cast to T&
 *
 *  these properties are tested here for each class defined in the namespace "buffered"
 */


TEST(BufferedFields, NdBitset){
    const std::array<size_t, 3> shape {4, 3, 2};
    const size_t nsetbit = 8;

    //  1. constructable via argument forwarding to T::T
    buffered::NdBitset<size_t, 3> direct(shape);
    const field::NdBitset<size_t, 3>& base_cref = direct;
    ASSERT_EQ(direct.m_format.m_nelement, 4 * 3 * 2);

    //  2. compatible with a call to any public method of T without cast to T&
    for (size_t i=0ul; i<nsetbit; ++i) direct.set(hashing::in_range(i, 0, direct.m_format.m_nelement));

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_shape, shape);
    for (size_t i=0ul; i<nsetbit; ++i)
        ASSERT_TRUE(direct_cpy.get(hashing::in_range(i, 0, direct_cpy.m_format.m_nelement)));
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    buffered::NdBitset<size_t, 3> base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_format.m_shape, shape);
    for (size_t i=0ul; i<nsetbit; ++i)
        ASSERT_TRUE(base_cref_cpy.get(hashing::in_range(i, 0, direct_cpy.m_format.m_nelement)));
    ASSERT_EQ(direct, base_cref_cpy);

    //  5. copy-assignable from const B<T>&
    direct.zero();
    ASSERT_TRUE(direct.is_zero());
    direct = direct_cpy;
    ASSERT_EQ(direct, direct_cpy);

    //  6. copy-assignable from const T&
    direct.zero();
    ASSERT_TRUE(direct.is_zero());
    direct = base_cref;
    ASSERT_EQ(direct, base_cref);

    //  7. assignable via any other method T::operator=
    direct.zero();
    ASSERT_EQ(direct.nsetbit(), 0ul);
    defs::inds setinds = {0, 1, 5, 7};
    direct = {0, 1, 5, 7};
    ASSERT_EQ(direct.nsetbit(), setinds.size());
    for (const auto& ibit : setinds) ASSERT_TRUE(direct.get(ibit));
}



TEST(BufferedFields, Numbers){
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
