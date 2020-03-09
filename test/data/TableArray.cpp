//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#include <src/data/Specification.h>
#include "src/data/Table.h"
#include "src/data/TableArray.h"
#include "src/data/BitfieldNew.h"

#if 0
TEST(TableArray, SingleTableCase) {
    Specification spec;
    const size_t nint = 4;
    const size_t nbit = 20;
    spec.add<size_t>(nint);
    spec.add<BitfieldNew>(nbit);
    const size_t nrow = 6;
    const size_t ntable = 1;

    TableArray<Table> table_array(ntable, spec, nrow);
    *table_array[0].view<size_t>(3) = 8;
    ASSERT_EQ(*table_array[0].view<size_t>(3), 8);
}


TEST(TableArray, EncodeDecode) {
    Specification spec;
    const size_t nint = 4;
    const size_t nbit = 20;
    spec.add<size_t>(nint);
    spec.add<BitfieldNew>(nbit);
    const size_t nrow = 6;
    const size_t ntable = 4;

    TableArray<Table> table_array(ntable, spec, nrow);
    *table_array[3].view<size_t>(3) = 8;
    ASSERT_EQ(*table_array[3].view<size_t>(3), 8);
}


#endif