//
// Created by Robert John Anderson on 2020-04-12.
//

#include <src/core/fermion/DecodedDeterminant.h>
#include <src/core/thread/Atomic.h>
#include "gtest/gtest.h"
#include "src/core/thread/PrivateStore.h"
#include "src/core/util/defs.h"

struct TestType1 {
    char i1, i2, i3;
    int i4;
    std::complex<float> z1;
    std::complex<double> z2, z3;
};

TEST(PrivateStore, AlignmentTest1) {
    ASSERT_EQ(sizeof(TestType1), 48);
    PrivateStore<TestType1> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    bool all_passed = true;
#pragma omp parallel default(none) shared(store, all_passed)
    {
        bool passed = true;
        /*
         * test whether every element is aligned to the cache line
         */
        if ((size_t) &store.get() % defs::cache_line_size) passed = false;
        as_atomic(all_passed) &= passed;
    }
    ASSERT_TRUE(all_passed);
}

struct TestType2 {
    char i1, i2, i3;
    std::complex<double> z1, z2;
    int i4, i5;
    std::complex<float> z3, z4;
    std::complex<double> z5, z6, z7;
    size_t i6, i7, i8;
};

TEST(PrivateStore, AlignmentTest2) {
    /*
     * stretches across multiple cacheline blocks
     */
    ASSERT_EQ(sizeof(TestType2), 136);
    PrivateStore<TestType2> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    bool all_passed = true;
#pragma omp parallel default(none) shared(store, all_passed)
    {
        bool passed = true;
        /*
         * test whether every element is aligned to the cache line
         */
        if ((size_t) &store.get() % defs::cache_line_size) passed = false;
        as_atomic(all_passed) &= passed;
    }
    ASSERT_TRUE(all_passed);
}

TEST(PrivateStore, ReduceLogical) {
    PrivateStore<bool> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = omp_get_thread_num() ? false : true;
    }
    ASSERT_EQ(omp_get_max_threads()==1, store.reduce_land());
    ASSERT_TRUE(store.reduce_lor());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = false;
    }
    ASSERT_FALSE(store.reduce_land());
    ASSERT_FALSE(store.reduce_lor());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = true;
    }
    ASSERT_TRUE(store.reduce_land());
    ASSERT_TRUE(store.reduce_lor());
}


TEST(PrivateStore, IntegralReduceSum) {
    PrivateStore<size_t> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = omp_get_thread_num();
    }
    ASSERT_EQ(store.reduce_sum(), (store.nthread() * (store.nthread() - 1)) / 2);
}

TEST(PrivateStore, IntegralReduceProd) {
    PrivateStore<size_t> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = omp_get_thread_num() + 1;
    }
    ASSERT_EQ(store.reduce_prod(), integer_utils::factorial(store.nthread()));
}

TEST(PrivateStore, IntegralReduceMax) {
    PrivateStore<size_t> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = (omp_get_thread_num() * 23) % 5;
    }
    size_t max = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = (i * 23) % 5;
        if (tmp > max) max = tmp;
    }
    ASSERT_EQ(store.reduce_max(), max);
}

TEST(PrivateStore, IntegralReduceMin) {
    PrivateStore<size_t> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = (omp_get_thread_num() * 23) % 5;
    }
    size_t min = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = (i * 23) % 5;
        if (tmp < min) min = tmp;
    }
    ASSERT_EQ(store.reduce_min(), min);
}

TEST(PrivateStore, RealReduceSum) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = 0.23 * omp_get_thread_num();
    }
    ASSERT_FLOAT_EQ(store.reduce_sum(), 0.23 * (store.nthread() * (store.nthread() - 1)) / 2);
}

TEST(PrivateStore, RealReduceProd) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = 0.23 * (omp_get_thread_num() + 1);
    }
    ASSERT_FLOAT_EQ(store.reduce_prod(),
                    std::pow(0.23, store.nthread()) * integer_utils::factorial(store.nthread()));
}

TEST(PrivateStore, RealReduceMax) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double((omp_get_thread_num() * 23) % 5) - 3;
    }
    double max = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = std::abs(double((i * 23) % 5) - 3);
        if (tmp > max) max = tmp;
    }
    ASSERT_FLOAT_EQ(store.reduce_max(), max);
}

TEST(PrivateStore, RealReduceMin) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double((omp_get_thread_num() * 23) % 5) - 3;
    }
    double min = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = std::abs(double((i * 23) % 5) - 3);
        if (tmp < min) min = tmp;
    }
    ASSERT_FLOAT_EQ(store.reduce_min(), min);
}

TEST(PrivateStore, ComplexReduceSum) {
    PrivateStore<std::complex<double>> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    const auto z = complex_utils::normal_from_sector(1, 5);
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double(omp_get_thread_num()) * z;
    }
    ASSERT_TRUE(consts::floats_equal(store.reduce_sum(),
                                     z * double((store.nthread() * (store.nthread() - 1)) / 2)));
}

TEST(PrivateStore, ComplexReduceProd) {
    PrivateStore<std::complex<double>> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    const auto z = complex_utils::normal_from_sector(1, 5);
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double(omp_get_thread_num()+1) * z;
    }
    ASSERT_TRUE(consts::floats_nearly_equal(
            store.reduce_prod(), std::pow(z, store.nthread()) *
            (double)integer_utils::factorial(store.nthread()), 1e-6));
}

TEST(PrivateStore, ComplexReduceMax) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    const auto z = complex_utils::normal_from_sector(1, 5);
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double((omp_get_thread_num() * 23) % 5) - 3;
    }
    double max = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = std::abs(double((i * 23) % 5) - 3);
        if (tmp > max) max = tmp;
    }
    ASSERT_FLOAT_EQ(store.reduce_max(), max);
}

TEST(PrivateStore, ComplexReduceMin) {
    PrivateStore<double> store;
    ASSERT_EQ(store.nthread(), omp_get_max_threads());
    const auto z = complex_utils::normal_from_sector(1, 5);
#pragma omp parallel default(none) shared(store)
    {
        store.get() = double((omp_get_thread_num() * 23) % 5) - 3;
    }
    double min = 0ul;
    for (size_t i = 0ul; i < store.nthread(); ++i) {
        auto tmp = std::abs(double((i * 23) % 5) - 3);
        if (tmp < min) min = tmp;
    }
    ASSERT_FLOAT_EQ(store.reduce_min(), min);
}
