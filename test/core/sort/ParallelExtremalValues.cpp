//
// Created by Robert John Anderson on 2020-08-02.
//

#include <src/core/fieldz/FieldsZ.h>
#include <src/core/fieldz/BufferedTableZ.h>
#include "gtest/gtest.h"
#include "src/core/sort/ParallelExtremalValues.h"
#include "src/core/sample/PRNG.h"


TEST(ParallelExtremalValues, Test) {

    typedef fieldsz::Number<double> field_t;
    struct TestRow: RowZ {
        field_t m_field;
        TestRow() : m_field(this){}
    };
    BufferedTableZ<TestRow> m_table("Test", {});
    PRNG prng(14, 1000);

    const size_t nrow_per_rank = 100;
    m_table.resize(nrow_per_rank);
    m_table.push_back(nrow_per_rank);

    auto &row = m_table.m_row;
    for (row.restart(); row.m_i < nrow_per_rank; row.step())
        row.m_field = prng.draw_float();

    /*
    const size_t nfind = 8;
    ParallelExtremalValues<field_t> pxv(m_table.m_row.m_field);
    pxv.reset(m_table);

    pxv.find(nfind);
     */

}