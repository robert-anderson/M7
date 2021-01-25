//
// Created by rja on 21/10/2020.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/Table.h"
#include "src/core/field/Fields.h"
#include "src/core/field/Elements.h"
#if 0
#include "src/core/data/BufferedField.h"



//template<size_t nind>
//struct ConfigurationField : Field {
//    NdFieldX<DeterminantField, 1> det;
//    NdFieldX<BosonOnvSpecifier, 1> perm;
//    ConfigurationField(DeterminantField&& det_field, BosonOnvSpecifier&& boson_field)
//};

struct FlagBase;

struct FlagField : BitsetSpecifier {
    TableX* m_table;
    std::vector<FlagBase*> m_flags;
    FlagField(TableX* table):BitsetSpecifier(0), m_table(table){}
    size_t add_flag(FlagBase* flag);
};

struct FlagBase {
    FlagField* m_spec;
    const size_t m_nbit;
    FlagBase(FlagField* field, size_t nbit):
    m_spec(field), m_nbit(nbit){}
    BitsetSpecifier::View::BitView operator()(const size_t& irow, const size_t& ibit){
        ASSERT(ibit<m_nbit);
        return BitsetSpecifier::View::BitView(
                BitsetSpecifier::View(*m_spec, m_spec->m_table->begin(irow)), ibit);
    }
};

size_t FlagBase::add_flag(FlagBase* flag){

}

//template<size_t nind>
//struct Flag : FlagBase {
//    Flag(FlagField* field):FlagBase(field){}
//};


struct MyFlags : FlagField {
    //Flag<1> deterministic;
    FlagBase deterministic;
    MyFlags(TableX* table): FlagField(table),
    deterministic(this, 3){}
};


using namespace fields;
struct TestTable : TableX {
    //FermionOnvSpecifier<0> field;
    Numbers<int, 1> ints;
    FermionOnv dets;
    NumberArray<double, 2> matrices;

    TestTable() :
            ints(this, "asdasdsadas", 3),
            dets(this, {10}, "asdasdsadas"),
            matrices(this, {10, 2}, "asdasdsadas")
            {}
};


TEST(Table_ND, Packing) {
    BufferedTable<TestTable> bt;
    //ASSERT_EQ(bt.m_row_dsize, 2);
    bt.print_column_details();
    bt.expand(10);
    bt.dets(0)[0] = true;
    bt.dets(0)[5] = true;
    std::cout << bt.dets(0).to_string() << std::endl;
}

//void f(const BufferedTable<TestTable> &bt) {
//    std::cout << bt.bits(0, 0).to_string() << std::endl;
//}



struct TestTable : TableX {
    fields::Onv config;
    TestTable():config(this, 4, 5, "Onv"){}
};

TEST(Table_ND, Test) {

    //BufferedField<FermionOnvSpecifier, 0> work_det({6}, {});
    elements::Onv work_config_buffer(6ul, 5ul);
    auto work_config = work_config_buffer();
    std::cout << work_config.to_string() << std::endl;
    }

#endif
//    BufferedTable<TestTable> bt;
//    bt.print_column_details();
//    bt.expand(10);
//    auto view = bt.config(0);
//    std::cout << view.to_string() << std::endl;


    //std::cout << typeid(fields::Numeric<int, 2>).name() << std::endl;
   // std::cout << typeid(fields::Numeric<int>).name() << std::endl;

//    BufferedTable<TestTable> bt;
//    bt.print_column_details();
//    bt.expand(10);
    //std::cout << bt.ints(0, 0) << std::endl;
//    f(bt);
    //BitsetSpecifier<0> field(&t, {}, 10, "asdasdsadas");
    //FermionOnvSpecifier<0> dets(&t, {}, 5, "asdasdsadas");
    //BosonOnvSpecifier<0> perms(&t, {}, 8, "zddsaas");
//    perms(0)[1] = 12;
//    perms(0)[4] = 11;
//
    //auto view = field(0);
    //std::cout << view.to_string() << std::endl;
//    std::cout << view.nboson() << std::endl;

//    std::vector<size_t> vec(6);
//    auto p = (size_t*)((char*)vec.data()+1);
//    p[0] = ~0ul;
//    for (auto i: vec) std::cout << i << std::endl;


//    Table t;
//    std::vector<double> vec(50, 9.99);
//    t.m_data = (char*)vec.data();
//    NumericArraySpecifier<double, 1, 2> matrices(&t, {6}, {2, 3}, "some matrices");
//    matrices(0, 0)(0, 0) = 0;
//    matrices(0, 0)(0, 1) = 1;
//    matrices(0, 0)(0, 2) = 2;
//    matrices(0, 0)(1, 0) = 3;
//    matrices(0, 0)(1, 1) = 4;
//    matrices(0, 0)(1, 2) = 5;
//    matrices(0, 2) = matrices(0, 0);
//    matrices.raw_view(0, 0).second;

    //for (auto item: vec) std::cout << item << std::endl;

