//
// Created by Robert John Anderson on 2020-04-12.
//

#include <src/core/fermion/DecodedDeterminant.h>
#include "gtest/gtest.h"
#include "src/core/thread/PrivateStore.h"
#include "src/defs.h"

struct alignas(defs::cache_line_size) TestType {
    char i1, i2;
    int i3;
    std::complex<float> z1;
    std::complex<double> z2, z3;
};

TEST(PrivateStore, Test) {
    const size_t nelement = 6;
    PrivateStore<TestType> store(nelement, TestType());
    bool all_passed = true;
#pragma omp parallel
    {
        bool passed = true;
        for (size_t ielement=0ul; ielement<nelement; ++ielement) {
            auto tmp = store.get(ielement);
            /*
             * test whether every element is aligned to the cache line
             */
            if ((size_t)&tmp % defs::cache_line_size) passed=false;
        }
#pragma omp atomic update
        all_passed &= passed;
    }
    ASSERT_TRUE(all_passed);
};