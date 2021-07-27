//
// Created by rja on 08/05/2021.
//

#include <src/core/table/BufferedFields.h>
#include "gtest/gtest.h"
#include "src/core/basis/CiSpaces.h"

TEST(CiSpaces, Test){
    typedef SingleFieldRow<fields::Mbf> row_t;
    const size_t nsite = 6;
    const size_t nelec = 6;
    BufferedTable<row_t> table("", {nsite});
    ci_gen::SpinSym gen(nsite, nelec, -1);
    gen(table.m_row, table.m_row.m_field);
    std::cout << table.to_string() << std::endl;
}