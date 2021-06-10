//
// Created by rja on 02/10/2020.
//

#include <src/core/table/BufferedFields.h>
#include "src/core/table//BufferedTable.h"
#include "gtest/gtest.h"

TEST(BufferedTable, Empty) {
    typedef SingleFieldRow<fields::Number<double>> row_t;
    BufferedTable<row_t> table("", {{}});
    ASSERT_EQ(table.m_nrow, 0);
    ASSERT_EQ(table.m_hwm, 0);
    auto& row = table.m_row;
    for (row.restart(); row.in_range(); row.step()){}
}

//
//
//struct TestTable1 : public Table_NEW {
//    NumericField<short, 1> shorts;
//    NumericField<float, 1> floats;
//    TestTable1(size_t n1, size_t n2):
//        shorts(this, "some shorts", n1),
//        floats(this, "some floats", n2)
//    {}
//};
//
//TEST(BufferedTable, DataIntegrityNumeric){
//    const size_t nshort = 7;
//    const size_t nfloat = 12;
//    BufferedTable<TestTable1> bt(nshort, nfloat);
//    const size_t nrow_ = 15;
//    bt.expand(nrow_);
//    for (size_t irow = 0ul; irow<nrow_; ++irow){
//        ASSERT_EQ(bt.push_back(), irow);
//        for (size_t ishort=0ul; ishort<nshort; ++ishort)
//            bt.shorts(irow, ishort) = ishort*nshort;
//        for (size_t ifloat=0ul; ifloat<nfloat; ++ifloat)
//            bt.floats(irow, ifloat) = ifloat*nfloat;
//    }
//    for (size_t irow = 0ul; irow<nrow_; ++irow){
//        for (size_t ishort=0ul; ishort<nshort; ++ishort)
//            ASSERT_EQ(ishort*nshort, bt.shorts(irow, ishort));
//        for (size_t ifloat=0ul; ifloat<nfloat; ++ifloat)
//            ASSERT_EQ(ifloat*nfloat, bt.floats(irow, ifloat));
//    }
//    ASSERT_EQ(bt.m_hwm, nrow_);
//    bt.clear();
//    ASSERT_EQ(bt.m_nrow, nrow_);
//    ASSERT_EQ(bt.m_hwm, 0);
//}
//
//struct TestTable2 : public Table_NEW {
//    NumericField<double, 1> doubles;
//    NumericField<unsigned char, 1> chars;
//    NumericVectorField<float, 1> float_vectors;
//    TestTable2(size_t ndouble, size_t nchar, size_t nvector, size_t nvector_item):
//    doubles(this, "some doubles", ndouble),
//    chars(this, "some chars", nchar),
//    float_vectors(this, nvector_item, "some vectors of floats", nvector)
//    {}
//};

//TEST(BufferedTable, DataIntegrityNumericVector) {
//    BufferedTable<TestTable2> bt;
//    bt.expand()
//    const size_t nshort = 7;
//    const size_t nfloat = 12;
//    BufferedTable<TestTable1> bt(nshort, nfloat);
//    const size_t nrow_ = 15;
//    bt.expand(nrow_);
//    for (size_t irow = 0ul; irow<nrow_; ++irow){
//}