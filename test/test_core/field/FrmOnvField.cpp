//
// Created by anderson on 09/02/2022.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"

TEST(FrmOnvField, SetFromInds) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    defs::inds setbits = {1, 90, nsite-1, nsite, 2*nsite-1};
    mbf = setbits;
    for (size_t ibit=0ul; ibit<mbf.m_nspinorb; ++ibit){
        bool is_set = std::find(setbits.cbegin(), setbits.cend(), ibit)!=setbits.cend();
        ASSERT_EQ(mbf.get(ibit), is_set);
    }
}

TEST(FrmOnvField, ClrSpinChannel) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), 0ul);
    mbf.put_spin_channel(0, true);
    ASSERT_EQ(mbf.nalpha(), nsite);
    ASSERT_EQ(mbf.nsetbit(), nsite);
    mbf.put_spin_channel(1, true);
    ASSERT_EQ(mbf.nalpha(), nsite);
    ASSERT_EQ(mbf.nsetbit(), 2*nsite);
    mbf.put_spin_channel(0, false);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), nsite);
    mbf.put_spin_channel(1, false);
    ASSERT_EQ(mbf.nalpha(), 0ul);
    ASSERT_EQ(mbf.nsetbit(), 0ul);
}

TEST(FrmOnvField, ForeachSetBit) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hashing::unique_in_range(0, 64, 0, mbf.m_nspinorb, true);
    mbf = setbits;

    auto it = setbits.cbegin();
    auto fn = [&it](size_t ibit){
        ASSERT_EQ(ibit, *(it++));
    };
    mbf.foreach_setbit(fn);
    // make sure all bits were iterated over
    ASSERT_EQ(it, setbits.cend());
}

TEST(FrmOnvField, ForeachSetBitPair) {
    const size_t nsite = 123;
    buffered::FrmOnv mbf(nsite);
    auto setbits = hashing::unique_in_range(0, 64, 0, mbf.m_nspinorb, true);
    mbf = setbits;

    auto it = setbits.cbegin();
    auto fn = [&it](size_t ibit, size_t jbit){
        std::cout << ibit << " " << jbit << std::endl;
        //ASSERT_EQ(ibit, *(it++));
    };
    mbf.foreach_setbit_pair(fn);
    // make sure all bits were iterated over
    //ASSERT_EQ(it, setbits.cend());
}