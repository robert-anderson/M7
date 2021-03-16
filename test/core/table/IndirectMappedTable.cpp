//
// Created by rja on 11/03/2021.
//

#include <src/core/table/MappedTable.h>
#include "gtest/gtest.h"

#if 0
TEST(IndirectMappedTable, Test) {
    struct IntRow : Row {
        fields::Number<size_t> m_number;
        IntRow(): m_number(this){}
        fields::Number<size_t> &key_field() {
            return m_number;
        };
    };

    MappedTable<IntRow> source({}, 10);
    //source.insert(4);

    //IndirectMappedTable
}
#endif