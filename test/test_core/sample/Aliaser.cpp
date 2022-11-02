//
// Created by Robert John Anderson on 2020-02-20.
//

#include "test_core/defs.h"
#include <M7_lib/defs.h>
#include <M7_lib/sample/Aliaser.h>

TEST(Aliaser, DistributionCheck) {
    PRNG prng(18, 1e4);

    v_t<prob_t> probs =
        {0, 0.648, 0.025, 0.035, 0.036, 0, 0.0648, 0.053, 0.0723, 0.1234};
    SingleAliaser aliaser(probs);
    uintv_t results(probs.size(), 0ul);

    const uint_t n_attempts = 20000000;
    for (uint_t i = 0ul; i < n_attempts; ++i) {
        results[aliaser.draw(prng)]++;
    }
    auto norm = std::accumulate(probs.begin(), probs.end(), 0.0);
    for (uint_t i = 0ul; i < probs.size(); ++i) {
        if (fptol::near_zero(probs[i])) {
            ASSERT_NEAR_ZERO(results[i]);
        }
        else {
            ASSERT_LT(std::abs(1.0 - (float(results[i]) / float(n_attempts)) / (probs[i] / norm)), 0.01);
        }
    }
}
