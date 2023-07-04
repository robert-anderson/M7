//
// Created by Robert J. Anderson on 1/7/22.
//

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>

#include "Hash.h"

#include <set>

using namespace hash;

digest_t hash::digest_in_range(digest_t v, digest_t lo, digest_t hi) {
    REQUIRE_GT(hi, lo, "upper hash value does not exceed lower value");
    v = fnv(v + 312194ul);
    v %= hi - lo;
    return v + lo;
}

digest_t hash::digest_in_range(const v_t<digest_t> &v, digest_t lo, digest_t hi) {
    REQUIRE_FALSE(v.empty(), "there must be as least one hashing value");
    auto out = digest_in_range(v[0], lo, hi);
    for (uint_t i = 1ul; i < v.size(); ++i) out = digest_in_range(out + 4321 * v[i], lo, hi);
    return out;
}