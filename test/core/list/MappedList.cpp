//
// Created by Robert John Anderson on 2020-04-01.
//

#include <src/core/table/NumericField.h>
#include "gtest/gtest.h"
#include "src/core/list/MappedList.h"

struct TestMappedList : public MappedList<NumericElement<size_t>> {
    NumericField<size_t> key, value;
    TestMappedList(size_t nbucket, size_t nsegment=1):
    MappedList("test mapped list", key, nbucket), key(this), value(this) {}
};

TEST(MappedList, DataIntegrity) {
    const size_t nrow = 36;
    TestMappedList list(10, 1);
    list.expand(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push(i * i);
        list.value(irow) = i;
    }
    auto irow = list.lookup_irow(27 * 27);
    ASSERT_NE(irow, ~0ul);
    ASSERT_EQ((size_t)list.value(irow), 27);
}


/*
struct MultipleHashingColumnsList : public MappedList<NumericElement<size_t>> {
    NumericField<size_t> key, value;
    TestMappedList(size_t nbucket, size_t nsegment=1):
            MappedList("test mapped list", key, nbucket), key(this), value(this) {}
};

TEST(MappedList, MultipleHashingColumns) {

}
*/