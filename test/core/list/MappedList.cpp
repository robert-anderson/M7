//
// Created by Robert John Anderson on 2020-04-01.
//

#include <src/core/table/NumericField.h>
#include "gtest/gtest.h"
#include "src/core/list/MappedList.h"

struct TestMappedList : public MappedList<NumericElement<size_t>> {
    NumericField<size_t> key, value;
    TestMappedList(size_t nbucket_per_segment, size_t nsegment=1):
    MappedList(key, nbucket_per_segment, nsegment),
    key(this), value(this) {}
};

TEST(MappedList, DataIntegrity) {
    const size_t nrow = 36;
    TestMappedList list(10, 1);
    list.expand(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push(i * i);
        list.value.element(irow) = i;
    }
    list.print();
    list.print_map();
    auto irow = list.lookup(27 * 27);
    ASSERT_NE(irow, ~0ul);
    ASSERT_EQ((size_t)list.value.element(irow), 27);
}

TEST(MappedList, ThreadSafety) {
    const size_t nrow = 10;
    TestMappedList list(nrow, 1);
    list.expand(nrow);

#pragma omp parallel for default(none) shared(list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto key = i*i;
            auto mutex = list.key_mutex(key);
            size_t irow = list.push(mutex, key);
            list.value.element(irow) = i;
        }
    }
    for (size_t i = 0; i < nrow; ++i) {
        auto key = i*i;
        auto irow = list.lookup(key);
        ASSERT_NE(irow, ~0ul);
        ASSERT_EQ((size_t)list.value.element(irow), i);
    }
}