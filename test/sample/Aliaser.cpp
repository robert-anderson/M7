//
// Created by Robert John Anderson on 2020-02-20.
//

#include <gtest/gtest.h>
#include <src/sample/PRNG.h>
#include <src/defs.h>
#include <src/sample/Aliaser.h>
#include <src/utils.h>

TEST(Aliaser, DistributionCheck){
    std::vector<defs::prob_t> probs = {2.3, 0.2, 0.5, 0.15, 0.15, 0.4, 1.2};
    Aliaser aliaser(probs);
    defs::inds results(probs.size(), 0ul);
    
    PRNG prng(0);
    const size_t n_attempts = 1000000;
    for (size_t i=0ul; i<n_attempts; ++i){
        results[aliaser.draw(prng)]++;
    }
    auto norm = std::accumulate(probs.begin(), probs.end(), 0.0);
    for (size_t i=0ul; i<probs.size(); ++i){
        ASSERT_LT(std::abs(1.0-(float(results[i])/float(n_attempts))/(probs[i]/norm)), 0.01);
    }
}