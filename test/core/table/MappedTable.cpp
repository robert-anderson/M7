//
// Created by rja on 08/10/2020.
//

#include <src/core/field/Elements.h>
#include "src/core/table/MappedTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/BufferedTable.h"
#include "gtest/gtest.h"

namespace mapped_table_test {
//    struct FloatTable : MappedTable<fields::Det> {
//        fields::Number<float> m_f;
//        FloatTable():
//    };

    struct DetTable : Table {
        fields::Det m_config;
        DetTable(size_t nsite) :
                m_config(this, nsite, "configuration") {}
    };

    struct DetMappedTable : MappedTable<DetTable, fields::Det> {
        DetMappedTable(size_t nsite):
        MappedTable<DetTable, fields::Det>(m_config, 10){}
    };
}

TEST(MappedTable, TEST) {
    const size_t nsite = 10;
    typedef BufferedTable<mapped_table_test::DetMappedTable> table_t;

    table_t bt("Mapped table test", nsite);
    bt.expand(10);
    elements::Det config(nsite);


//    config.m_fonv[2] = 1;
//    config.m_bonv[2] = 5;
//    std::cout << config.to_string() << std::endl;
//    bt.insert(config);
//    config.m_bonv[4] = 6;
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
//    std::cout << bt.field_details() << std::endl;
//    //bt.insert(bitset);
//}