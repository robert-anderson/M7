//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#include <M7_lib/defs.h>
#include <M7_lib/sample/Aliaser.h>

TEST(Aliaser, DistributionCheck) {
    PRNG prng(18, 1e4);

    std::vector<defs::prob_t> probs =
        {0, 0.648, 0.025, 0.035, 0.036, 0, 0.0648, 0.053, 0.0723, 0.1234};
    Aliaser aliaser(probs);
    defs::inds results(probs.size(), 0ul);

    const size_t n_attempts = 20000000;
    for (size_t i = 0ul; i < n_attempts; ++i) {
        results[aliaser.draw(prng)]++;
    }
    auto norm = std::accumulate(probs.begin(), probs.end(), 0.0);
    for (size_t i = 0ul; i < probs.size(); ++i) {
        if (consts::nearly_zero(probs[i])) ASSERT_EQ(results[i], 0);
        else ASSERT_LT(std::abs(1.0 - (float(results[i]) / float(n_attempts)) / (probs[i] / norm)), 0.01);
    }
}