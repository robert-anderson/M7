//
// Created by rja on 21/03/2021.
//

#include "gtest/gtest.h"
#include "M7_lib/enumerator/Enumerator.h"

TEST(Enumerator, CombinationsDistinct){
    const size_t n=5, r=2;
    enums::CombinationsDistinct e(n, r);
    std::vector<defs::inds> correct = {
            {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4},
            {2, 3}, {2, 4},
            {3, 4}
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}

TEST(Enumerator, CombinationsDistinctCache){
    const size_t n=5, r=2;
    enums::Cache<enums::CombinationsDistinct> e(n, r);
    std::vector<defs::inds> correct = {
            {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 2}, {1, 3}, {1, 4},
            {2, 3}, {2, 4},
            {3, 4}
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}

TEST(Enumerator, CombinationsWithRepetition){
    const size_t n=5, r=2;
    enums::CombinationsWithRepetition e(n, r);
    std::vector<defs::inds> correct = {
            {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 1}, {1, 2}, {1, 3}, {1, 4},
            {2, 2}, {2, 3}, {2, 4},
            {3, 3}, {3, 4},
            {4, 4}
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}

TEST(Enumerator, CombinationsWithRepetitionCache){
    const size_t n=5, r=2;
    enums::Cache<enums::CombinationsWithRepetition> e(n, r);
    std::vector<defs::inds> correct = {
            {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 1}, {1, 2}, {1, 3}, {1, 4},
            {2, 2}, {2, 3}, {2, 4},
            {3, 3}, {3, 4},
            {4, 4}
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}

TEST(Enumerator, PermutationsWithRepetition){
    const size_t n=5, r=2;
    enums::PermutationsWithRepetition e(n, r);
    std::vector<defs::inds> correct = {
            {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
            {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4},
            {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 4},
            {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4},
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}

TEST(Enumerator, PermutationsWithRepetitionCache){
    const size_t n=5, r=2;
    enums::Cache<enums::PermutationsWithRepetition> e(n, r);
    std::vector<defs::inds> correct = {
            {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
            {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
            {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4},
            {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 4},
            {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4},
    };
    auto c = correct.cbegin();
    while (e.next()){
        for (size_t i=0; i<r; ++i) ASSERT_EQ(e[i], (*c)[i]);
        ++c;
    };
}
