//
// Created by rja on 08/10/2020.
//

#include <src/core/field/Elements.h>
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "gtest/gtest.h"

#if 0
struct TestTable : MappedTable<fields::Configuration> {
    fields::Configuration m_config;
    TestTable(size_t nsite, size_t nmode):
    MappedTable<fields::Configuration>(m_config, 100),
            m_config(this, nsite, nmode, "configuration"){}
};

TEST(MappedTable, TEST){
    const size_t nsite = 10;
    const size_t nmode = 8;
    BufferedTable<TestTable> bt(nsite, nmode);
    bt.expand(10);
    elements::Configuration config(nsite, nmode);
    auto config_view = config();
    config_view.m_det[2] = 1;
    std::cout << config_view.to_string() << std::endl;
    bt.insert(config_view);
    bt.erase(bt[config_view]);
    auto lookup = bt[config_view];
    std::cout << bool(lookup) << std::endl;

}

#endif

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
//    std::cout << bt.field_details() << std::endl;
//    //bt.insert(bitset);
//}