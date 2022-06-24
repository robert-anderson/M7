//
// Created by Robert John Anderson on 2020-02-21.
//

#include <gtest/gtest.h>
#include "M7_lib/sample/PRNG.h"

TEST(PRNG, MeanCheck) {
    PRNG prng(0, 1000);
    const uint_t n=10000000;
    double tot = 0;
    for (uint_t i=0; i<n; ++i){
        tot+=prng.draw_float();
    }
    ASSERT_EQ(std::round(1000*2*tot/n), 1000);
}

TEST(PRNG, StochasticRound) {
    PRNG prng(0, 1000);
    const uint_t n=10000000;
    double v = 123.34;
    double rounding_magnitude = 1.3;
    double tot = 0;
    for (uint_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_LT(std::abs(tot/n-v), 1e-4);
}

TEST(PRNG, NegativeStochasticRound) {
    PRNG prng(0, 1000);
    const uint_t n=10000000;
    double v = -123.34;
    double rounding_magnitude = 1.3;
    double tot = 0;
    for (uint_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_LT(std::abs(tot/n-v), 1e-4);
}

TEST(PRNG, ComplexStochasticRound) {
    PRNG prng(0, 1000);
    const uint_t n=10000000;
    std::complex<double> v = {-123.34, 54.2};
    double rounding_magnitude = 0.4;
    std::complex<double> tot = 0;
    for (uint_t i=0; i<n; ++i){
        tot+=prng.stochastic_round(v, rounding_magnitude);
    }
    ASSERT_LT(std::abs(tot/(double)n-v), 1e-4);
}

TEST(PRNG, StochasticThreshold) {
    PRNG prng(0, 1000);
    double rounding_magnitude = 2.7;
    ASSERT_EQ(prng.stochastic_threshold(123.34, rounding_magnitude), 123.34);
    prob_t accept_prob;
    prng.stochastic_threshold(2.5, rounding_magnitude, accept_prob);
    ASSERT_FLOAT_EQ(accept_prob, 2.5 / 2.7);

    const uint_t n=20000000;
    double v = -0.14;
    double tot = 0;
    for (uint_t i=0; i<n; ++i){
        tot+=prng.stochastic_threshold(v, rounding_magnitude);
    }
    ASSERT_LT(std::abs(tot/n-v), 1e-4);
}

TEST(PRNG, ModularBase){
    PRNG prng(0, 1000);
    const uint32_t modular_base = 6;
    std::vector<uint_t> frequencies(modular_base, 0ul);
    const uint_t n=10000000;
    for (uint_t i=0; i<n; ++i){
        frequencies[prng.draw_uint(modular_base)]++;
    }
    ASSERT_TRUE(std::all_of(frequencies.cbegin(), frequencies.cend(), [](const uint_t &i){return i>0;}));
}

TEST(PRNG, InRange) {
    PRNG prng(0, 100000);
    const uint32_t min = 5, max = 10;
    const uint_t n = 100000000;
    uint_t tot = 0ul;
    for (uint_t i = 0; i < n; ++i) {
        tot += prng.draw_uint(min, max);
    }
    ASSERT_LT(std::abs(tot / double(n) - (max+min-1)/2.0), 1e-4);
}