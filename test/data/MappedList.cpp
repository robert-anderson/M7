//
// Created by Robert John Anderson on 2020-02-23.
//

#include <gtest/gtest.h>
#include "src/data/MappedList.h"

#if 0
TEST(MappedList, DataIntegrity) {
    Specification specification;
    specification.add<size_t>(2);
    size_t nrow = 36;

    MappedList<size_t> list(specification, nrow, 0);

    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push(i * i);
        list.view<size_t>(irow)[1] = i;
    }
    list.print();
    auto irow = list.lookup(27 * 27);
    ASSERT_EQ(list.view<size_t>(irow)[1], 27);
}

TEST(MappedList, ThreadSafety) {
    Specification specification;
    specification.add<size_t>(2);
    size_t nrow = 36;

    MappedList<size_t> list(specification, nrow, 0);

#pragma omp parallel for default(none) shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto mutex = list.key_mutex(i * i);
            size_t irow = list.push(mutex, i * i);
            list.view<size_t>(irow)[1] = i;
        }
    }
    for (size_t i = 0; i < nrow; ++i) {
        auto irow = list.lookup(i*i);
        ASSERT_EQ(list.view<size_t>(irow)[1], i);
    }
}
#endif