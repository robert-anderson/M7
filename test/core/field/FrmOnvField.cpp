//
// Created by anderson on 09/02/2022.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"

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

TEST(FrmOnvField, ClearSpinChannel) {
//    const size_t nsite = 123;
//    buffered::FrmOnv mbf(nsite);
//
//    const size_t i
//    mbf.set(ibegin, iend);
//

//    for (size_t ibit=0ul; ibit<mbf.m_nspinorb; ++ibit) mbf.set(ibit);
//    ASSERT_EQ(mbf.nalpha(), nsite);
//    ASSERT_EQ(mbf.nsetbit(), 2*nsite);
}