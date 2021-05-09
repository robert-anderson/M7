//
// Created by Robert John Anderson on 2020-08-02.
//

#include <src/core/field/Fields.h>
#include <src/core/table/BufferedTable.h>
#include <src/core/table/BufferedFields.h>
#include <src/core/sort/LocalExtremalValues.h>
#include "gtest/gtest.h"
#include "src/core/sort/ParallelExtremalValues.h"
#include "src/core/sample/PRNG.h"


TEST(ParallelExtremalValues, Test) {

    typedef SingleFieldRow<fields::Number<size_t>> row_t;
    BufferedTable<row_t> table("Test", {{}});
    const size_t nrow_per_rank = 40;

    auto get_value = [](size_t irank, size_t ielement){
        return hashing::in_range((irank+3)*(ielement+9), 0, nrow_per_rank*mpi::nrank());
    };

    table.push_back(nrow_per_rank);
    auto row = table.m_row;
    for (row.restart(); row.in_range(); row.step())
        row.m_field = get_value(mpi::irank(), row.m_i);
    ASSERT_EQ(table.m_hwm, nrow_per_rank);

    const size_t nfind = 10;
    LocalExtremalValues<row_t, size_t, 0ul> lxv(table.m_row, table.m_row.m_field, 1, 1);
    lxv.find(nfind);
    for (size_t i=0ul; i<lxv.m_xinds.nfound(); ++i) {
        row.jump(lxv.m_xinds[i]);
        std::cout << row.m_field << std::endl;
    }

    /*
    auto row_cmp = row;
    auto cmp_fn = [&](const size_t& irow, const size_t& irow_cmp){
        row.jump(irow);
        row_cmp.jump(irow_cmp);
        return row.m_field < row_cmp.m_field;
    };

    ExtremalIndices xv(cmp_fn);

    const size_t nfind = 10;
    xv.reset(table.m_hwm);
    xv.find(nfind);
    for (size_t i=0ul; i<xv.nfound(); ++i){
        row.jump(xv[i]);
        std::cout << row.m_field << std::endl;
    }
    */
std::cout << "" << std::endl;
    if (mpi::i_am_root()){
        defs::inds all_values;
        all_values.reserve(nrow_per_rank*mpi::nrank());
        for (size_t irank=0ul; irank<mpi::nrank(); ++irank)
            for (size_t ielement=0ul; ielement<nrow_per_rank; ++ielement)
                all_values.emplace_back(get_value(irank, ielement));
        std::sort(all_values.begin(), all_values.end());

        for (size_t i=0ul; i<nfind; ++i)
            std::cout << all_values[i] << std::endl;
    }

    /*
    const size_t nfind = 8;
    ParallelExtremalValues<field_t> pxv(m_table.m_row.m_field);
    pxv.reset(m_table);

    pxv.find(nfind);
     */

}