//
// Created by rja on 03/03/2021.
//

#include "src/core/table/BufferedFields.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/sort/QuickSorter.h"
#include "src/core/wavefunction/SpawnTable.h"
#include "src/core/wavefunction/WalkerTable.h"
#include "gtest/gtest.h"

TEST(QuickSorter, Test){
    const size_t nsite = 5;
    const size_t nrow = 10;

    BufferedTable<SpawnTableRow> table("test", {{nsite, false}});
    table.push_back(nrow);
    auto row = table.m_row;

    /*
     * defining a stock of arbitrary MBFs
     */
    buffered::FrmOnv onv168(nsite);
    onv168 = {1, 6, 8};
    buffered::FrmOnv onv239(nsite);
    onv239 = {2, 3, 9};
    buffered::FrmOnv onv345(nsite);
    onv345 = {3, 4, 5};
    buffered::FrmOnv onv468(nsite);
    onv468 = {4, 6, 8};

    row.restart();
    row.m_dst_mbf = onv168;
    row.step(); row.m_dst_mbf = onv239;
    row.step(); row.m_dst_mbf = onv345;
    row.step(); row.m_dst_mbf = onv345;
    row.step(); row.m_dst_mbf = onv239;
    row.step(); row.m_dst_mbf = onv345;
    row.step(); row.m_dst_mbf = onv239;
    row.step(); row.m_dst_mbf = onv345;
    row.step(); row.m_dst_mbf = onv168;
    row.step(); row.m_dst_mbf = onv468;

    std::vector<fields::FrmOnv*> correct_order = {&onv239, &onv345, &onv168, &onv468};
    std::vector<size_t> correct_counts = {3, 4, 2, 1};

    auto row1 = table.m_row;
    auto row2 = table.m_row;
    auto comp_fn = [&](const size_t &irow1, const size_t &irow2){
        row1.jump(irow1);
        row2.jump(irow2);
        return row1.m_dst_mbf < row2.m_dst_mbf;
    };

    /*
     * sorting in ascending lexical order
     */
    QuickSorter qs(comp_fn);
    qs.preserve_sort(table);

    /*
     * check that table is sorted into proper blocks
     */
    auto correct = correct_order.cbegin();
    for (size_t i=0ul; i<table.m_hwm; ++i){
        row.jump(qs.m_inds[i]);
        while (correct!=correct_order.cend() && row.m_dst_mbf != **correct) correct++;
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
