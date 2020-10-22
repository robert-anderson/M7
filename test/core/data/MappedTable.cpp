//
// Created by rja on 08/10/2020.
//
//
//#include "src/core/data/MappedTable.h"
//#include "src/core/data/BufferedTable.h"
//#include "src/core/data/NumericField.h"
//#include "src/core/data/NumericVectorField.h"
//#include "src/core/data/BitsetField.h"
//#include "gtest/gtest.h"
//
//struct TestTable : public MappedTable<1> {
//    BitsetField<0> bitset;
//    NumericField<size_t, 0> index;
//    NumericVectorField<size_t, 0> vector;
//    TestTable(): MappedTable(100, &bitset),
//    bitset(this, 10, "some bitset"),
//    index(this, "some index"),
//    vector(this, 12, "some vector")
//    {}
//};
//
//TEST(MappedTable, SingleFieldInsert){
//    BufferedTable<TestTable> bt;
//    bt.expand(100);
//    Bitset<0> bitset(10);
//    std::cout << bitset.to_string() << std::endl;
//    bt.push_back();
//    bt.bitset(0);// = bitset;
//    std::cout << bt.bitset(0).to_string() << std::endl;
//    std::cout << bt.field_details() << std::endl;
//    //bt.insert(bitset);
//}