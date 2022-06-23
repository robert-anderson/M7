//
// Created by rja on 23/06/22.
//

#ifndef M7_PROBABILITY_H
#define M7_PROBABILITY_H

#include <numeric>
#include "M7_lib/defs.h"

/**
 * Contains methods for managing simple unit-normalized vectors of probabilities.
 * For structures which can be probed in constant time, see sample/Aliaser
 */
namespace prob {

    static void normalize(std::vector <defs::prob_t> &v, defs::prob_t norm = 1.0) {
        auto tot = std::accumulate(v.begin(), v.end(), 0.0);
        auto fac = norm / tot;
        for (auto &i: v) i *= fac;
    }

    static void rectify(std::vector <defs::prob_t> &v, defs::prob_t min) {
        for (auto &prob: v) if (prob < min) prob = min;
        normalize(v);
    }
}

#endif //M7_PROBABILITY_H
