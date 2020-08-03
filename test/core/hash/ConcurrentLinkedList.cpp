//
// Created by Robert John Anderson on 2020-08-01.
//

#include "gtest/gtest.h"
#include <omp.h>
#include "src/core/hash/ConcurrentLinkedList.h"

TEST(LinkedList, Test) {
    ConcurrentLinkedList<size_t> list;
    ASSERT_TRUE(list.is_empty());
    const size_t n = 100;
#pragma omp parallel for
    for (size_t i=0ul; i<n; ++i){
        list.append(i*i);
    }
    ASSERT_EQ(list.size(), n);

    auto second = list.m_next->m_next;
    list.delete_after(&list);
    ASSERT_EQ(second, list.m_next);
}