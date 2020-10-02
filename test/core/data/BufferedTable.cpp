//
// Created by rja on 02/10/2020.
//

#include "src/core/data/BufferedTable.h"
#include "src/core/data/NumericField.h"
#include "gtest/gtest.h"


struct TestTable1 : public Table {
    NumericField<short, 1> shorts;
    NumericField<float, 1> floats;
    TestTable1(size_t n1, size_t n2):
        shorts(this, n1),
        floats(this, n2)
    {}
};

TEST(BufferedTable, DataIntegrityNumeric){
    const size_t nshort = 7;
    const size_t nfloat = 12;
    BufferedTable<TestTable1> bt(nshort, nfloat);
    const size_t nrow = 15;
    bt.expand(nrow);
    for (size_t irow = 0ul; irow<nrow; ++irow){
        ASSERT_EQ(bt.push_back(), irow);
        for (size_t ishort=0ul; ishort<nshort; ++ishort)
            bt.shorts(irow, ishort) = ishort*nshort;
        for (size_t ifloat=0ul; ifloat<nfloat; ++ifloat)
            bt.floats(irow, ifloat) = ifloat*nfloat;
    }
    for (size_t irow = 0ul; irow<nrow; ++irow){
        for (size_t ishort=0ul; ishort<nshort; ++ishort)
            ASSERT_EQ(ishort*nshort, bt.shorts(irow, ishort));
        for (size_t ifloat=0ul; ifloat<nfloat; ++ifloat)
            ASSERT_EQ(ifloat*nfloat, bt.floats(irow, ifloat));
    }
    ASSERT_EQ(bt.m_hwm, nrow);
    bt.clear();
    ASSERT_EQ(bt.m_nrow, nrow);
    ASSERT_EQ(bt.m_hwm, 0);
}