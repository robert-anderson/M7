//
// Created by rja on 08/10/2020.
//

#include "gtest/gtest.h"
#include "TableTest.h"

TEST(MappedTable, TEST) {
//    const size_t nsite = 10;
//    typedef BufferedTable<table_test::DetMappedTable> table_t;
//
//    table_t bt("Mapped table test", nsite);
//    bt.expand(10);
//    elements::Det config(nsite);



//    config.m_frm[2] = 1;
//    config.m_bos[2] = 5;
//    std::cout << config.to_string() << std::endl;
//    bt.insert(config);
//    config.m_bos[4] = 6;
//    bt.insert(config);
////    bt.erase(bt[config]);
////    auto lookup = bt[config];
////    std::cout << bool(lookup) << std::endl;
//    bt.print_map();

}


//#include "src/core/data/BufferedTable.h"
//#include "src/core/data/NumericField.h"
//#include "src/core/data/NumericVectorField.h"
//#include "src/core/data/BitsetField.h"
//
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
//    std::cout << bt.column_details() << std::endl;
//    //bt.insert(bitset);
//}