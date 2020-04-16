//
// Created by Robert John Anderson on 2020-02-21.
//

#include <gtest/gtest.h>
#include "src/core/sample/PRNG.h"

TEST(PRNG, MeanCheck) {
    PRNG prng(0, 100);
    const size_t n=10000000;
    double tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.draw_float();
    }
    ASSERT_EQ((int)(1000*2*tot/n), 1000);
}

TEST(PRNG, StochasticRound) {
    PRNG prng(0, 100);
    const size_t n=50000000;
    double v = 123.34;
    double rounding_magnitude = 1.3;
    double tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_TRUE(std::abs(tot/n-v)<1e-4);
}

TEST(PRNG, NegativeStochasticRound) {
    PRNG prng(0, 100);
    const size_t n=10000000;
    double v = -123.34;
    double rounding_magnitude = 3.7;
    double tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_TRUE(std::abs(tot/n-v)<1e-4);
}

TEST(PRNG, ComplexStochasticRound) {
    PRNG prng(0, 100);
    const size_t n=50000000;
    std::complex<double> v = {-123.34, 54.2};
    double rounding_magnitude = 0.4;
    std::complex<double> tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_TRUE(std::abs(tot/(double)n-v)<1e-4);
}

TEST(PRNG, StochasticThreshold) {
    PRNG prng(0, 100);
    double rounding_magnitude = 3.7;
    ASSERT_EQ(prng.stochastic_threshold(123.34, rounding_magnitude), 123.34);

    const size_t n=10000000;
    double v = -0.34;
    double tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.stochastic_threshold(v, rounding_magnitude);
    }
    ASSERT_TRUE(std::abs(tot/n-v)<1e-4);
}
