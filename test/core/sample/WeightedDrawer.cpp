//
// Created by RJA on 20/11/2020.
//

#include "gtest/gtest.h"
#include "src/core/sample/WeightedDrawer.h"

TEST(WeightedDrawer, Test) {
    PRNG prng(14, 1000);
    WeightedDrawer drawer(4, prng);
    drawer.set(0.2, 0.4, 0.05);
    ASSERT_FLOAT_EQ(drawer.prob(3), 0.35);
    drawer.set(0.2, 0.4, 0.1, 0.3);
    defs::inds freqs(4, 0ul);
    const size_t ndraw = 100000000;
    for (size_t idraw=0ul; idraw<ndraw; ++idraw){
        freqs[drawer.draw()]++;
    }
    ASSERT_TRUE(consts::floats_nearly_equal(0.2/(freqs[0]/(double)ndraw), 1.0, 1e-3));
    ASSERT_TRUE(consts::floats_nearly_equal(0.4/(freqs[1]/(double)ndraw), 1.0, 1e-3));
    ASSERT_TRUE(consts::floats_nearly_equal(0.1/(freqs[2]/(double)ndraw), 1.0, 1e-3));
    ASSERT_TRUE(consts::floats_nearly_equal(0.3/(freqs[3]/(double)ndraw), 1.0, 1e-3));
}