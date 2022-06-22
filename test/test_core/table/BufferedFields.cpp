//
// Created by Robert J. Anderson on 09/06/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/basis/BasisData.h"


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
    const auto elements = hash::unique_in_range(0, nelement, 0, hi, true);

    //  1. constructable via argument forwarding to T::T
    buffered::Numbers<hash::digest_t, 3> direct(shape);
    ASSERT_EQ(direct.m_format.m_nelement, nelement);

    //  2. compatible with a call to any public method of T without cast to T&
    for (size_t ielement = 0ul; ielement<nelement; ++ielement) direct[ielement] = elements[ielement];

    //  3. copy-constructable from const B<T>&
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_format.m_shape, shape);
    ASSERT_EQ(direct.to_vector(), elements);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::Numbers<hash::digest_t, 3>& base_cref = direct;
    buffered::Numbers<hash::digest_t, 3> base_cref_cpy(base_cref);
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
    const auto setbits = hash::unique_in_range(0, nsetbit, 0, nbit, true);

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
    auto setbits = hash::unique_in_range(0, nsetbit, 0, nbit, true);

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


TEST(BufferedFields, FrmOnv){
    const size_t nsite = 6;
    defs::inds_t frm_inds = {0, 2, 4, 5, 7, 9};
    const sys::frm::Basis basis(nsite);

    //  1. constructable via argument forwarding to T::T
    buffered::FrmOnv direct(basis);
    ASSERT_EQ(direct.m_basis, basis);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.nalpha(), 0ul);

    //  3. copy-constructable from const B<T>&
    direct = frm_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_basis, direct.m_basis);
    ASSERT_EQ(direct_cpy.nsetbit(), frm_inds.size());
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::FrmOnv& base_cref = direct;
    ASSERT_EQ(base_cref.nsetbit(), frm_inds.size());
    buffered::FrmOnv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_basis, direct.m_basis);
    ASSERT_EQ(base_cref_cpy.nsetbit(), frm_inds.size());
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
    direct = frm_inds;
    ASSERT_EQ(direct, direct_cpy);
}

TEST(BufferedFields, BosOnv){
    const size_t nmode = 5;
    const size_t bos_occ_cutoff = 10;
    defs::inds_t zero_inds(nmode, 0);
    defs::inds_t bos_inds = {3, 1, 2, 4, 9};
    const sys::bos::Basis basis(nmode, bos_occ_cutoff);

    //  1. constructable via argument forwarding to T::T
    buffered::BosOnv direct(nmode);
    ASSERT_EQ(direct.m_basis.m_nmode, nmode);
    ASSERT_EQ(direct.to_vector<size_t>(), zero_inds);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_nelement, nmode);

    //  3. copy-constructable from const B<T>&
    direct = bos_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_basis, direct.m_basis);
    ASSERT_EQ(direct_cpy.to_vector<size_t>(), bos_inds);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::BosOnv& base_cref = direct;
    buffered::BosOnv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_basis, direct.m_basis);
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
    direct = bos_inds;
    ASSERT_EQ(direct, direct_cpy);
}

TEST(BufferedFields, FrmBosOnv){
    defs::inds_t frm_inds = {0, 2, 4, 5, 7};
    defs::inds_t bos_inds = {3, 1, 2, 4, 9};
    const sys::frm::Basis frm_basis(4);
    const sys::bos::Basis bos_basis(5, defs::max_bos_occ);

    //  1. constructable via argument forwarding to T::T
    buffered::FrmBosOnv direct(frm_basis, bos_basis);
    ASSERT_EQ(direct.m_frm.m_basis, frm_basis);
    ASSERT_EQ(direct.m_bos.m_basis, bos_basis);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_bos.m_nelement, bos_basis.m_nmode);
    ASSERT_EQ(direct.m_frm.nalpha(), 0ul);

    //  3. copy-constructable from const B<T>&
    direct.m_frm = frm_inds;
    direct.m_bos = bos_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_frm.m_basis, direct.m_frm.m_basis);
    ASSERT_EQ(direct_cpy.m_bos.m_basis, direct.m_bos.m_basis);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::FrmBosOnv& base_cref = direct;
    buffered::FrmBosOnv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_frm.m_basis, direct.m_frm.m_basis);
    ASSERT_EQ(base_cref_cpy.m_bos.m_basis, direct.m_bos.m_basis);
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

TEST(BufferedFields, FrmXonv){
    const size_t nsite = 5;
    defs::inds_t ket_inds = {0, 2, 4, 5, 9};
    defs::inds_t bra_inds = {1, 2, 3, 7};
    const sys::frm::Basis basis(nsite);

    //  1. constructable via argument forwarding to T::T
    buffered::FrmXonv direct(basis);
    ASSERT_EQ(direct.m_ket.m_basis, basis);
    ASSERT_EQ(direct.m_bra.m_basis, basis);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_ket.nalpha(), 0ul);
    ASSERT_EQ(direct.m_bra.nalpha(), 0ul);

    //  3. copy-constructable from const B<T>&
    direct.m_ket = ket_inds;
    direct.m_bra = bra_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_ket.m_basis, direct.m_ket.m_basis);
    ASSERT_EQ(direct_cpy.m_bra.m_basis, direct.m_ket.m_basis);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::FrmXonv& base_cref = direct;
    buffered::FrmXonv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_ket.m_basis, direct.m_ket.m_basis);
    ASSERT_EQ(base_cref_cpy.m_bra.m_basis, direct.m_bra.m_basis);
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
    direct = {ket_inds, bra_inds};
    ASSERT_EQ(direct, direct_cpy);
}

TEST(BufferedFields, BosXonv){
    const size_t nmode = 5;
    defs::inds_t ket_inds = {3, 1, 2, 4, 9};
    defs::inds_t bra_inds = {2, 4, 9, 3, 1};
    const sys::bos::Basis basis(nmode);

    //  1. constructable via argument forwarding to T::T
    buffered::BosXonv direct(nmode);
    ASSERT_EQ(direct.m_ket.m_basis, basis);
    ASSERT_EQ(direct.m_bra.m_basis, basis);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_ket.m_nelement, nmode);
    ASSERT_EQ(direct.m_bra.m_nelement, nmode);

    //  3. copy-constructable from const B<T>&
    direct.m_ket = ket_inds;
    direct.m_bra = bra_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_ket.m_basis, direct.m_ket.m_basis);
    ASSERT_EQ(direct_cpy.m_bra.m_basis, direct.m_bra.m_basis);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::BosXonv& base_cref = direct;
    buffered::BosXonv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_ket.m_basis, direct.m_ket.m_basis);
    ASSERT_EQ(base_cref_cpy.m_bra.m_basis, direct.m_bra.m_basis);
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
    direct = {ket_inds, bra_inds};
    ASSERT_EQ(direct, direct_cpy);
}

TEST(BufferedFields, FrmBosXonv){
    defs::inds_t ket_frm_inds = {0, 2, 4, 5, 7};
    defs::inds_t ket_bos_inds = {3, 1, 2, 4, 9};
    defs::inds_t bra_frm_inds = {1, 2, 3, 7};
    defs::inds_t bra_bos_inds = {2, 4, 9, 3, 1};
    const sys::frm::Basis frm_basis(4);
    const sys::bos::Basis bos_basis(5);

    //  1. constructable via argument forwarding to T::T
    buffered::FrmBosXonv direct(frm_basis, bos_basis);
    ASSERT_EQ(direct.m_ket.m_frm.m_basis, frm_basis);
    ASSERT_EQ(direct.m_ket.m_bos.m_basis, bos_basis);
    ASSERT_EQ(direct.m_bra.m_frm.m_basis, frm_basis);
    ASSERT_EQ(direct.m_bra.m_bos.m_basis, bos_basis);

    //  2. compatible with a call to any public method of T without cast to T&
    ASSERT_EQ(direct.m_ket.m_bos.m_nelement, bos_basis.m_nmode);
    ASSERT_EQ(direct.m_ket.m_frm.nalpha(), 0ul);
    ASSERT_EQ(direct.m_bra.m_bos.m_nelement, bos_basis.m_nmode);
    ASSERT_EQ(direct.m_bra.m_frm.nalpha(), 0ul);

    //  3. copy-constructable from const B<T>&
    direct.m_ket.m_frm = ket_frm_inds;
    direct.m_ket.m_bos = ket_bos_inds;
    direct.m_bra.m_frm = bra_frm_inds;
    direct.m_bra.m_bos = bra_bos_inds;
    auto direct_cpy = direct;
    ASSERT_EQ(direct_cpy.m_ket.m_frm.m_basis, direct.m_ket.m_frm.m_basis);
    ASSERT_EQ(direct_cpy.m_ket.m_bos.m_basis, direct.m_ket.m_bos.m_basis);
    ASSERT_EQ(direct_cpy.m_bra.m_frm.m_basis, direct.m_bra.m_frm.m_basis);
    ASSERT_EQ(direct_cpy.m_bra.m_bos.m_basis, direct.m_bra.m_bos.m_basis);
    ASSERT_EQ(direct, direct_cpy);

    //  4. copy-constructable from const T&
    const field::FrmBosXonv& base_cref = direct;
    buffered::FrmBosXonv base_cref_cpy(base_cref);
    ASSERT_EQ(base_cref_cpy.m_ket.m_frm.m_basis, direct.m_ket.m_frm.m_basis);
    ASSERT_EQ(base_cref_cpy.m_ket.m_bos.m_basis, direct.m_ket.m_bos.m_basis);
    ASSERT_EQ(base_cref_cpy.m_bra.m_frm.m_basis, direct.m_bra.m_frm.m_basis);
    ASSERT_EQ(base_cref_cpy.m_bra.m_bos.m_basis, direct.m_bra.m_bos.m_basis);
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
    direct = {{ket_frm_inds, ket_bos_inds}, {bra_frm_inds, bra_bos_inds}};
    ASSERT_EQ(direct, direct_cpy);
}
