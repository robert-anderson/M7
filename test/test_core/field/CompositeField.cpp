//
// Created by anderson on 1/25/22.
//

#include "gtest/gtest.h"
#include "src/core/field/CompositeField.h"
#include "src/core/table/BufferedFields.h"


struct DetPerm : CompositeField<field::FrmOnv, field::BosOnv> {
    typedef CompositeField<field::FrmOnv, field::BosOnv> base_t;
    field::FrmOnv m_frm;
    field::BosOnv m_bos;
    DetPerm(Row* row, BasisData bd, std::string prefix):
        base_t(m_frm, m_bos),
        m_frm(row, bd.m_nsite, prefix+" fermions"),
        m_bos(row, bd.m_nmode, prefix+" bosons"){}

    DetPerm(const DetPerm& other): base_t(m_frm, m_bos), m_frm(other.m_frm), m_bos(other.m_bos){}

    DetPerm& operator=(const DetPerm& other) {
        m_frm = other.m_frm;
        m_bos = other.m_bos;
        return *this;
    }
};

struct OuterDetPerm : CompositeField<DetPerm, DetPerm> {
    typedef CompositeField<DetPerm, DetPerm> base_t;
    DetPerm m_ket;
    DetPerm m_bra;
    OuterDetPerm(Row* row, BasisData bd): base_t(m_ket, m_bra),
        m_ket(row, bd, "ket"), m_bra(row, bd, "bra"){}

    OuterDetPerm(const OuterDetPerm& other): base_t(m_ket, m_bra), m_ket(other.m_ket), m_bra(other.m_bra){}

    OuterDetPerm& operator=(const OuterDetPerm& other) {
        m_ket = other.m_ket;
        m_bra = other.m_bra;
        return *this;
    }

};

struct CompFieldTestRow : Row {
    field::FrmOnv m_frm_onv;
    OuterDetPerm m_dmbas;
    CompFieldTestRow() :
            m_frm_onv(this, 6, "alfa"),
            m_dmbas(this, {8, 4}){}
};

struct DetPermTestRow : Row {
    DetPerm m_dp;
    DetPermTestRow() : m_dp(this, {6, 4}, "www"){}
};

//struct NestedCompFieldTestRow : Row {
//    field::FrmOnv m_frm_onv;
//    typedef CompositeField<field::FrmOnv, field::FrmOnv> frm_onv_pair_t;
//    frm_onv_pair_t m_comp_field;
//    CompositeField<frm_onv_pair_t, field::BosOnv> m_nested_field;
//
//    NestedCompFieldTestRow() :
//            m_frm_onv({this, 6, "alfa"}),
//            m_comp_field({this, 6, "bravo"}, {this, 8, "charlie"}),
//            m_nested_field({{this, 10, "delta"}, {this, 5, "echo"}},
//                           {this, 9, "foxtrot"}){}
//};


TEST(CompositeField, Test) {

//    buffered::BosOnv b(5);
//    b = {1, 5, 6, 9, 19};
//    buffered::BosOnv c(5);
//    c = static_cast<const field::BosOnv&>(b);
    buffered::FrmBosOnv b({6, 4});
    //buffered::FrmBosOnv c({6, 4});
    field::FrmBosOnv& bref(b);
    b = {{1, 5, 7}, {1, 5, 6, 9}};
    buffered::FrmBosOnv c(bref);
    std::cout << b.to_string() << std::endl;
    std::cout << c.to_string() << std::endl;
    c = b;
    std::cout << c.to_string() << std::endl;



//    DetPermTestRow row;
//    ASSERT_EQ(&row.m_dp.m_frm, &row.m_dp.get<0>());
}

#if 0

TEST(CompositeField, CopyMoveSemantics) {

    CompFieldTestRow row;
    ASSERT_EQ(&row.m_dmbas.m_ket.m_frm, &row.m_dmbas.m_ket.get<0>());

    ASSERT_EQ(row.nfield(), 5ul);
    for (size_t i=0ul; i<row.nfield(); ++i) ASSERT_EQ(row.m_fields[i]->m_row, &row);
    ASSERT_EQ(row.m_child, nullptr);

    auto row_copy = row;
    for (size_t i=0ul; i<row_copy.nfield(); ++i) ASSERT_EQ(row_copy.m_fields[i]->m_row, &row_copy);

    BufferedTable<CompFieldTestRow> bt("test table", {{}});
    bt.resize(10);
    bt.push_back(4);

    auto r1 = bt.m_row;
    r1.restart();
    r1.m_dmbas.m_ket.m_frm = {0, 1, 2, 3};
    r1.m_dmbas.m_bra.m_bos = {1, 2, 3, 7};
    auto r2 = bt.m_row;
    r2.jump(2);
    r2.m_dmbas = r1.m_dmbas;
    ASSERT_TRUE(r1.begin());
    ASSERT_TRUE(r2.begin());

    ASSERT_EQ(r1.m_frm_onv.m_row, &r1);
    ASSERT_EQ(r1.m_dmbas.m_ket.m_frm.m_row, &r1);
    ASSERT_EQ(r2.m_frm_onv.m_row, &r2);
    ASSERT_EQ(r2.m_dmbas.m_ket.m_frm.m_row, &r2);

    ASSERT_EQ(&r2.m_dmbas.m_ket.m_frm, &r2.m_dmbas.m_ket.get<0>());

    std::cout << (r2.m_frm_onv==r1.m_frm_onv) << std::endl;
    std::cout << (r2.m_dmbas==r1.m_dmbas) << std::endl;
    std::cout << (r2.m_dmbas!=r1.m_dmbas) << std::endl;
}

//TEST(CompositeField, NestedCopyMoveSemantics) {
//
//    //NestedCompFieldTestRow row;
//
//    BufferedTable<NestedCompFieldTestRow> bt("test table", {{}});
//    std::cout << bt.m_row.m_fields.size() << std::endl;
//    bt.resize(10);
//    bt.push_back(4);
//
//    auto row = bt.m_row;
//    row.restart();
//    row.m_nested_field.get<0>().get<0>() = {1, 3, 5};
//    std::cout << bt.to_string() << std::endl;
//
//}
#endif