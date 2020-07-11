//
// Created by Robert John Anderson on 2020-04-01.
//

#include <src/core/table/NumericField.h>
#include "gtest/gtest.h"
#include "src/core/list/MappedList.h"

struct TestMappedList : public MappedList<NumericElement<size_t>> {
    NumericField<size_t> key, value;
    TestMappedList(size_t nbucket, size_t nsegment=1):
    MappedList(key, nbucket), key(this), value(this) {}
};

TEST(MappedList, DataIntegrity) {
    const size_t nrow = 36;
    TestMappedList list(10, 1);
    list.expand(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push(i * i);
        list.value(irow) = i;
    }
    auto irow = list.lookup(27 * 27);
    ASSERT_NE(irow, ~0ul);
    ASSERT_EQ((size_t)list.value(irow), 27);
}

TEST(MappedList, ThreadSafety) {
    const size_t nrow = 3600;
    TestMappedList list(nrow/10, 1);
    list.expand(nrow);

#pragma omp parallel for default(none) shared(list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto key = i*i;
            auto mutex = list.key_mutex(key);
            size_t irow = list.push(mutex, key);
            list.value(irow) = i;
        }
    }
    for (size_t i = 0; i < nrow; ++i) {
        auto key = i*i;
        auto irow = list.lookup(key);
        ASSERT_NE(irow, ~0ul);
        ASSERT_EQ((size_t)list.value(irow), i);
    }
}

TEST(MappedList, ThreadSerialization) {
    const size_t nrow = 1e7;
    TestMappedList list(100000, 1);
    list.expand(nrow);

    auto do_something = [](size_t i){while (i%213){i = 3*i+1;}};

#pragma omp parallel for default(none) shared(list, do_something)
    for (size_t i = 0; i < nrow; ++i) {
        do_something(i);
        auto key = i*i;
        auto mutex = list.key_mutex(key);
        size_t irow = list.push(mutex, key);
        list.value(irow) = i;
    }
}