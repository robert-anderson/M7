//
// Created by Robert J. Anderson on 1/7/22.
//

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>

#include "Hash.h"

#include <set>

using namespace hash;

digest_t hash::in_range(digest_t v, digest_t lo, digest_t hi) {
    REQUIRE_GT(hi, lo, "upper hash value does not exceed lower value");
    v = fnv(v + 312194ul);
    v %= hi - lo;
    return v + lo;
}

digest_t hash::in_range(const v_t<digest_t> &v, digest_t lo, digest_t hi) {
    REQUIRE_FALSE(v.empty(), "there must be as least one hashing value");
    auto out = in_range(v[0], lo, hi);
    for (uint_t i = 1ul; i < v.size(); ++i) out = in_range(out + 4321 * v[i], lo, hi);
    return out;
}

v_t<digest_t> hash::in_range(const v_t<digest_t> &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted) {
    v_t<digest_t> out;
    out.reserve(ngen);
    auto vatt = v;
    vatt.push_back(0);
    while (out.size() != ngen) {
        auto r = in_range(vatt, lo, hi);
        out.push_back(r);
        ++vatt.back();
    }
    if (sorted) std::sort(out.begin(), out.end());
    return out;
}

v_t<digest_t> hash::in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted) {
    return in_range(v_t<digest_t>{v}, ngen, lo, hi, sorted);
}

v_t<digest_t> hash::unique_digest_in_range(const v_t<digest_t> &v, uint_t ngen, digest_t lo, digest_t hi, bool sorted) {
    REQUIRE_LE(lo + ngen, hi, "number of unique values can't exceed the range of allowed values");
    v_t<digest_t> out;
    out.reserve(ngen);
    std::set<digest_t> set;
    auto vatt = v;
    vatt.push_back(0);
    while (set.size() != ngen) {
        auto r = in_range(vatt, lo, hi);
        while (set.find(r) != set.end()) if (++r == hi) r = lo;
        out.push_back(r);
        set.insert(r);
        ++vatt.back();
    }
    if (sorted) std::sort(out.begin(), out.end());
    return out;
}

v_t<digest_t> hash::unique_digest_in_range(digest_t v, uint_t ngen, digest_t lo, digest_t hi, bool sorted) {
    return unique_digest_in_range(v_t<digest_t>{v}, ngen, lo, hi, sorted);
}
