//
// Created by rja on 09/06/2021.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"


/*
 * Explicitly implemented buffered fields B<T> of T (where T is any class derived from FieldBase or CompositeField) in
 * the namespace "buffered" must be:
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

TEST(BufferedFields, Numbers){
    const std::array<size_t, 3> shape {4, 3, 2};
    const size_t nelement = 4*3*2;
    const size_t hi = nelement*4;
    const auto elements = hashing::unique_in_range(0, nelement, 0, hi, true);

    //  1. constructable via argument forwarding to T::T
    buffered::Numbers<defs::hash_t, 3> direct(shape);
    ASSERT_EQ(direct.m_format.m_nelement, nelement);

    //  2. compatible with a call to any public method of T without cast to T&
    for (size_t ielement = 0ul; ielement<nelement; ++ielement) direct[ielement] = elements[ielement];

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_shape, shape);
    ASSERT_EQ(direct.to_vector(), elements);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::Numbers<defs::hash_t, 3>& base_cref = direct;
    buffered::Numbers<defs::hash_t, 3> base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_format.m_shape, shape);
    ASSERT_EQ(base_cref_cpy.to_vector(), elements);
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
    ASSERT_TRUE(direct.is_zero());
    direct = elements;
    ASSERT_EQ(direct.to_vector(), elements);
}


TEST(BufferedFields, Number){
    typedef std::complex<double> T;
    const T element = {1.98, -4.3};

    //  1. constructable via argument forwarding to T::T
    buffered::Number<T> direct;
    ASSERT_EQ(direct.m_format.m_nelement, 1ul);

    //  2. compatible with a call to any public method of T without cast to T&
    direct = element;

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_nelement, 1ul);
    ASSERT_EQ(direct.to_vector()[0], element);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::Number<T>& base_cref = direct;
    buffered::Number<T> base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_format.m_nelement, 1ul);
    ASSERT_EQ(base_cref_cpy.to_vector()[0], element);
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
    ASSERT_TRUE(direct.is_zero());
    direct = element;
    ASSERT_EQ(direct, element);
}

TEST(BufferedFields, NdBitset){
    const std::array<size_t, 3> shape {4, 3, 2};
    const size_t nbit = 4*3*2;
    const size_t nsetbit = 8;
    const auto setbits = hashing::unique_in_range(0, nsetbit, 0, nbit, true);


    //  1. constructable via argument forwarding to T::T
    buffered::NdBitset<size_t, 3> direct(shape);
    ASSERT_EQ(direct.m_format.m_nelement, nbit);

    //  2. compatible with a call to any public method of T without cast to T&
    for (auto& ibit: setbits) direct.set(ibit);

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_shape, shape);
    for (auto& ibit: setbits) ASSERT_TRUE(direct.get(ibit));
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::NdBitset<size_t, 3>& base_cref = direct;
    buffered::NdBitset<size_t, 3> base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_format.m_shape, shape);
    for (auto& ibit: setbits) ASSERT_TRUE(base_cref_cpy.get(ibit));
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
    direct = setbits;
    ASSERT_EQ(direct.nsetbit(), setbits.size());
    for (const auto& ibit : setbits) ASSERT_TRUE(direct.get(ibit));
}

TEST(BufferedFields, Bitset){
    const size_t nbit = 20;
    const size_t nsetbit = 8;
    auto setbits = hashing::unique_in_range(0, nsetbit, 0, nbit, true);

    //  1. constructable via argument forwarding to T::T
    buffered::Bitset<size_t> direct(nbit);
    ASSERT_EQ(direct.m_format.m_nelement, nbit);

    //  2. compatible with a call to any public method of T without cast to T&
    for (auto& ibit: setbits) direct.set(ibit);

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_nelement, nbit);
    for (auto& ibit: setbits) ASSERT_TRUE(direct.get(ibit));
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::Bitset<size_t>& base_cref = direct;
    buffered::Bitset<size_t> base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_format.m_nelement, nbit);
    for (auto& ibit: setbits) ASSERT_TRUE(direct.get(ibit));
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
    direct = setbits;
    ASSERT_EQ(direct.nsetbit(), setbits.size());
    for (const auto& ibit : setbits) ASSERT_TRUE(direct.get(ibit));
}


TEST(BufferedFields, FrmBosOnv){
    const BasisDims bd = {4, 5};
    defs::inds frm_inds = {0, 2, 4, 5, 7};
    defs::inds bos_inds = {3, 1, 2, 4, 9};

    //  1. constructable via argument forwarding to T::T
    buffered::FrmBosOnv direct(bd);
    ASSERT_EQ(direct.m_frm.m_nsite, bd.m_nsite);
    ASSERT_EQ(direct.m_bos.m_nmode, bd.m_nmode);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_bos.m_nelement, bd.m_nmode);
    ASSERT_EQ(direct.m_frm.nalpha(), 0ul);

    //  3. copy-constructable from const B<T>&
    direct.m_frm = frm_inds;
    direct.m_bos = bos_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_frm.m_nsite, direct.m_frm.m_nsite);
    ASSERT_EQ(direct_cpy.m_bos.m_nmode, direct.m_bos.m_nmode);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::FrmBosOnv& base_cref = direct;
    buffered::FrmBosOnv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_frm.m_nsite, direct.m_frm.m_nsite);
    ASSERT_EQ(base_cref_cpy.m_bos.m_nmode, direct.m_bos.m_nmode);
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
    direct = {frm_inds, bos_inds};
    ASSERT_EQ(direct, direct_cpy);
}
