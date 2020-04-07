//
// Created by Robert John Anderson on 2020-02-21.
//

#include <gtest/gtest.h>
#include "src/core/sample/PRNG.h"

TEST(PRNG, MeanCheck) {
    PRNG prng(0, 100);
    const size_t n=10000000;
    float tot = 0;
    for (size_t i=0; i<n; ++i){
        tot+=prng.draw_float();
    }
    ASSERT_EQ((int)(1000*2*tot/n), 1000);
}