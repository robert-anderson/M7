//
// Created by Robert John Anderson on 2020-02-09.
//


#include <gtest/gtest.h>
#include "src/data/Table.h"

TEST(Table, AllToAllV) {
    Specification spec;
    spec.create<int>(3);
    spec.create<float>(2);
    spec.create<BitfieldNew>(18);
    Table table(spec, 2, 3);
    *table.view<int>(1, 1, 2) = 1233456;
    table.bitfield_view(1, 0, 0).set(defs::inds{0, 2, 3, 4, 17});
    table.print();
}