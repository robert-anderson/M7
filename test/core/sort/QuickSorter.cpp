//
// Created by rja on 03/03/2021.
//

#include "src/core/table/BufferedFields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/sort/QuickSorter.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/dynamics/WalkerTable.h"
#include "gtest/gtest.h"

#ifndef ENABLE_BOSONS
TEST(QuickSorter, Test){
    const size_t nsite = 5;
    const size_t nrow = 10;

    BufferedTable<SpawnTableRow> table("test", {{nsite}});
    table.push_back(nrow);
    auto row = table.m_row;

    /*
     * defining a stock of arbitrary ONVs
     */
    buffered::Onv<> onv168(nsite);
    onv168 = {1, 6, 8};
    buffered::Onv<> onv239(nsite);
    onv239 = {2, 3, 9};
    buffered::Onv<> onv345(nsite);
    onv345 = {3, 4, 5};
    buffered::Onv<> onv468(nsite);
    onv468 = {4, 6, 8};

    row.restart();
    row.m_dst_onv = onv168;
    row.step(); row.m_dst_onv = onv239;
    row.step(); row.m_dst_onv = onv345;
    row.step(); row.m_dst_onv = onv345;
    row.step(); row.m_dst_onv = onv239;
    row.step(); row.m_dst_onv = onv345;
    row.step(); row.m_dst_onv = onv239;
    row.step(); row.m_dst_onv = onv345;
    row.step(); row.m_dst_onv = onv168;
    row.step(); row.m_dst_onv = onv468;

    std::vector<fields::Onv<>*> correct_order = {&onv239, &onv345, &onv168, &onv468};
    std::vector<size_t> correct_counts = {3, 4, 2, 1};

    std::cout <<
    table.to_string()
    << std::endl;


    //std::function<bool(const size_t &, const size_t &)> f = ;
    auto row1 = table.m_row;
    auto row2 = table.m_row;
    auto comp_fn = [&](const size_t &irow1, const size_t &irow2){
        row1.jump(irow1);
        row2.jump(irow2);
        return row1.m_dst_onv <= row2.m_dst_onv;
    };

    /*
     * sorting in ascending lexical order
     */
    Quicksorter qs(comp_fn);
    qs.sort(table);

    /*
     * check that table is sorted into proper blocks
     */
    auto correct = correct_order.cbegin();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        row.jump(qs.m_inds[i]);
        while (correct!=correct_order.cend() && row.m_dst_onv != **correct) correct++;
    }
    /*
     * correct ordering iterator should not be exhausted
     */
    ASSERT_FALSE(correct==correct_order.cend());
    /*
     * should be on last element of correct ordering so incrementation exhausts the iterator
     */
    ASSERT_TRUE(++correct==correct_order.cend());

    std::cout <<
              table.to_string(&qs.m_inds)
              << std::endl;
}
#endif