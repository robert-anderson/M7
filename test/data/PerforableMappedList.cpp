//
// Created by Robert John Anderson on 2020-02-23.
//

#include <gtest/gtest.h>
#include "src/data/PerforableMappedList.h"

#if 0
TEST(PerforableMappedList, Removal) {
    Specification specification;
    specification.add<size_t>(2);
    size_t nrow = 36;

    PerforableMappedList<size_t> list(specification, nrow, 0);

#pragma omp parallel for default(none) shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto mutex = list.key_mutex(i * i);
            size_t irow = list.push(mutex, i * i);
            list.view<size_t>(irow)[1] = i;
        }
    }
    list.remove(6*6);
    ASSERT_EQ(list.lookup(6*6), ~0ul);
}

TEST(PerforableMappedList, RemovalAndReuse) {
    Specification specification;
    specification.add<size_t>(2);
    size_t nrow = 10;

    PerforableMappedList<size_t> list(specification, nrow, 0);

#pragma omp parallel for default(none) shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto mutex = list.key_mutex(i * i);
            size_t irow = list.push(mutex, i * i);
            list.view<size_t>(irow)[1] = i;
        }
    }
    auto irow_available = list.remove(6*6);
    list.synchronize();
    ASSERT_EQ(list.push(nrow * nrow), irow_available);
}
#endif