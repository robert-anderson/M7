//
// Created by Robert John Anderson on 2020-08-02.
//

#include "gtest/gtest.h"
#include "src/core/sort/ParallelExtremalValues.h"
#include "src/core/field/Elements.h"
#include "src/core/sample/PRNG.h"


TEST(ParallelExtremalValues, Test) {

    typedef fields::Number<double> field_t;
    BufferedSingleFieldTable<field_t> m_table;
    PRNG prng(14, 1000);

    const size_t nrow_per_rank = 100;
    const size_t nfind = 8;
    m_table.resize(nrow_per_rank);
    m_table.push_back(nrow_per_rank);

    for (size_t i = 0ul; i < nrow_per_rank; ++i) m_table.m_field(i) = prng.draw_float();

    ParallelExtremalValues<field_t> pxv(m_table.m_field);
    pxv.reset(m_table);

    pxv.find(nfind);

}