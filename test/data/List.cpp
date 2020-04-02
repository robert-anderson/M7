//
// Created by Robert John Anderson on 2020-02-20.
//

#if 0
#include <gtest/gtest.h>
#include "src/data/Field.h"
#include "src/data/ListNew.h"

TEST(List, ThreadSafety) {
    struct TestList : public ListNew {
        Field<size_t> counter;
        TestList() :
                ListNew(),
                counter(this)
                {}
    };

    size_t nrow = 3600;
    TestList list;
    list.extend(nrow);

#pragma omp parallel for default(none), shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push();
        list.counter(irow) = irow;
    }

    for (auto irow{0ul}; irow < list.nrow(); ++irow) {
        ASSERT_EQ(list.counter(irow), irow);
    }
}

#endif