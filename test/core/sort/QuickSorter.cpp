//
// Created by rja on 03/03/2021.
//

#include "src/core/table/BufferedFields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/sort/LambdaQuickSorter.h"
#include "src/core/wavefunction/SpawnTable.h"
#include "src/core/wavefunction/WalkerTable.h"
#include "gtest/gtest.h"


TEST(QuickSorter, TableSort){
    const size_t nrow = 18;
    const size_t nchar = 9;

    typedef SingleFieldRow<field::String> row_t;

    BufferedTable<row_t> table("test", {{nchar}});
    table.push_back(nrow);
    auto& row = table.m_row;

    /*
     * fill table with random words
     */
    row.restart(); row.m_field = "alpha";
    row.step(); row.m_field = "beetle";
    row.step(); row.m_field = "catapult";
    row.step(); row.m_field = "alpaca";
    row.step(); row.m_field = "alpha";
    row.step(); row.m_field = "catapult";
    row.step(); row.m_field = "beetle";
    row.step(); row.m_field = "beetle";
    row.step(); row.m_field = "alpaca";
    row.step(); row.m_field = "alpha";
    row.step(); row.m_field = "catamaran";
    row.step(); row.m_field = "catapult";
    row.step(); row.m_field = "alpaca";
    row.step(); row.m_field = "catamaran";
    row.step(); row.m_field = "alpaca";
    row.step(); row.m_field = "alpaca";
    row.step(); row.m_field = "catapult";
    row.step(); row.m_field = "beetle";

    std::vector<std::string> correct_order = {"alpaca", "alpha", "beetle", "catamaran", "catapult"};
    std::vector<size_t> correct_counts = {5, 4, 3, 2, 4};

    auto row1 = table.m_row;
    auto row2 = table.m_row;
    auto comp_fn = [&](const size_t &irow1, const size_t &irow2){
        row1.jump(irow1);
        row2.jump(irow2);
        return row1.m_field < row2.m_field;
    };

    /*
     * sorting in ascending lexical order without physically reordering rows
     */
    LambdaQuickSorter qs(comp_fn);
    qs.preserve_sort(table);

    /*
     * check that table is sorted into proper blocks
     */
    auto correct = correct_order.cbegin();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        row.jump(qs.m_inds[i]);
        while (correct!=correct_order.cend() && row.m_field != *correct) correct++;
    }
    /*
     * correct ordering iterator should not be exhausted
     */
    ASSERT_FALSE(correct==correct_order.cend());
    /*
     * should be on last element of correct ordering so incrementation exhausts the iterator
     */
    ASSERT_TRUE(++correct==correct_order.cend());

    /*
     * sorting in ascending lexical order, this time physically reordering rows
     */
    qs.reorder_sort(table);
    correct = correct_order.cbegin();
    row.restart();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        while (correct!=correct_order.cend() && row.m_field != *correct) correct++;
        row.step();
    }
    /*
     * correct ordering iterator should not be exhausted
     */
    ASSERT_FALSE(correct==correct_order.cend());
    /*
     * should be on last element of correct ordering so incrementation exhausts the iterator
     */
    ASSERT_TRUE(++correct==correct_order.cend());

}