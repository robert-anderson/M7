//
// Created by Robert J. Anderson on 03/03/2021.
//

#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/table/BufferedTable.h"
#include "M7_lib/sort/LambdaQuickSorter.h"
#include "M7_lib/sort/QuickSorter.h"
#include "gtest/gtest.h"

namespace quick_sorter_test {

    typedef SingleFieldRow<field::String> row_t;
    typedef BufferedTable<row_t> bt_t;
    static void setup_table(Table<row_t>& table){
        const std::vector<std::string> words = {
                "alpha", "beetle", "catapult", "alpaca", "alpha", "catapult",
                "beetle", "beetle", "alpaca", "alpha", "catamaran", "catapult",
                "alpaca", "catamaran", "alpaca", "alpaca", "catapult", "beetle"
        };
        auto nrow = words.size();
        table.push_back(nrow);
        auto& row = table.m_row;
        /*
         * fill table with random words
         */
        row.restart();
        for (auto& word: words){
            row.m_field = word;
            row.step();
        }
    }

    static std::vector<std::string> correct_order() {
        return {"alpaca", "alpha", "beetle", "catamaran", "catapult"};
    }

    static std::vector<size_t> correct_counts() {
        return {5, 4, 3, 2, 4};
    }

    struct AscOrderFn {
        row_t m_row1;
        row_t m_row2;
        AscOrderFn(const row_t& row): m_row1(row), m_row2(row){}
        bool operator()(const size_t &irow1, const size_t &irow2){
            m_row1.jump(irow1);
            m_row2.jump(irow2);
            return m_row1.m_field < m_row2.m_field;
        };
    };
}

/*
 * test that the sorter works when the comparison behaviour is specified through statically dispatched code
 *  (should be faster than the Lambda equivalent)
 */
TEST(QuickSorter, FunctorSort){
    using namespace quick_sorter_test;
    const size_t nchar = 9;

    bt_t table("test", {{nchar}});
    setup_table(table);

    AscOrderFn comp_fn(table.m_row);
    /*
     * sorting in ascending lexical order without physically reordering rows
     */
    LambdaQuickSorter2 qs(comp_fn);
    qs.preserve_sort(table);

    auto& row1 = comp_fn.m_row1;
    /*
     * check that table is sorted into proper blocks
     */
    auto correct_order = quick_sorter_test::correct_order();
    auto correct = correct_order.cbegin();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        row1.jump(qs.m_inds[i]);
        while (correct!=correct_order.cend() && row1.m_field != *correct) correct++;
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
    row1.restart();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        while (correct!=correct_order.cend() && row1.m_field != *correct) correct++;
        row1.step();
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


/*
 * test that the sorter also works when the functor specified is actually a wrapper for dynamically dispatched code
 */
TEST(QuickSorter, LambdaSort){
    using namespace quick_sorter_test;
    const size_t nchar = 9;

    bt_t table("test", {{nchar}});
    setup_table(table);

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
    LambdaQuickSorter2 qs(comp_fn);
    qs.preserve_sort(table);

    /*
     * check that table is sorted into proper blocks
     */
    auto correct_order = quick_sorter_test::correct_order();
    auto correct = correct_order.cbegin();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        row1.jump(qs.m_inds[i]);
        while (correct!=correct_order.cend() && row1.m_field != *correct) correct++;
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
    row1.restart();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        while (correct!=correct_order.cend() && row1.m_field != *correct) correct++;
        row1.step();
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