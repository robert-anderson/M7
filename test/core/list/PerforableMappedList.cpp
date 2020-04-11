//
// Created by Robert John Anderson on 2020-02-23.
//

#include <gtest/gtest.h>
#include <src/core/list/PerforableMappedList.h>

struct TestPerforableMappedList : public PerforableMappedList<NumericElement<size_t>> {
    NumericField<size_t> key, value;
    TestPerforableMappedList(size_t nbucket):
    PerforableMappedList(key, nbucket), key(this), value(this) {}
};

TEST(PerforableMappedList, Removal) {
    const size_t nrow = 3600;
    TestPerforableMappedList list(nrow/10);
    list.expand(nrow);

#pragma omp parallel for default(none) shared(list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto mutex = list.key_mutex(i * i);
            size_t irow = list.push(mutex, i * i);
            list.value(irow) = i;
        }
    }
    list.remove(6*6);
    ASSERT_EQ(list.lookup(6*6), ~0ul);
}

TEST(PerforableMappedList, RemovalAndReuse) {
    const size_t nrow = 3600;
    TestPerforableMappedList list(nrow/10);
    list.expand(nrow);

#pragma omp parallel for default(none) shared(list)
    for (size_t i = 0; i < nrow; ++i) {
        {
            //enclose within scope to ensure destruction of mutex
            auto mutex = list.key_mutex(i * i);
            size_t irow = list.push(mutex, i * i);
            list.value(irow) = i;
        }
    }
    auto irow_available = list.remove(6*6);
    list.synchronize();
    ASSERT_EQ(list.push(nrow * nrow), irow_available);
}
