//
// Created by rja on 08/10/2020.
//

#include "src/core/data/MappedTable.h"
#include "src/core/data/BufferedTable.h"
#include "src/core/data/NumericField.h"
#include "src/core/data/BitsetField.h"
#include "gtest/gtest.h"

struct TestTable : public MappedTable<1> {
    BitsetField<0> bitset;
    NumericField<size_t, 0> index;
    TestTable(): MappedTable(100, &bitset), bitset(this, 6), index(this){}
};

TEST(MappedTable, SingleFieldInsert){
    //BufferedTable<TestTable> bt;
    //bt.expand(100);
    Bitset<0> bitset(10);
    std::cout << bitset.to_string() << std::endl;
    //bt.insert();
}