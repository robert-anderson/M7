//
// Created by Robert John Anderson on 2020-03-31.
//

#include <src/core/table/NumericField.h>
#include "gtest/gtest.h"
#include "src/core/list/List.h"

struct TestList : public List{
    NumericField<int> counter;
    TestList(size_t nsegment=1):
    List(nsegment),
    counter(this)
    {}
};


TEST(List, ThreadSafety) {
    size_t nrow = 3600;
    TestList list;
    list.expand(nrow);

#pragma omp parallel for default(none), shared(nrow, list)
    for (size_t i = 0; i < nrow; ++i) {
        size_t irow = list.push();
        list.counter.element(irow) = irow;
    }

    for (auto irow{0ul}; irow < list.nrow_per_segment(); ++irow) {
        ASSERT_EQ(list.counter.element(irow), irow);
    }
}