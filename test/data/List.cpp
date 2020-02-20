//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#include "src/data/List.h"


TEST(List, ThreadSafety) {
    Specification specification;
    specification.add<size_t>(1);
    size_t nrow = 3600;

    List list(specification, nrow);

#pragma omp parallel for
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push();
        list.view<size_t>(irow)[0] = irow;
    }

    for (auto irow{0ul}; irow < list.nrow(); ++irow) {
        auto v = list.view<size_t>(irow);
        ASSERT_EQ(*v, irow);
    }
}