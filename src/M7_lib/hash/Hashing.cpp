//
// Created by Robert J. Anderson on 1/7/22.
//

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>

#include "Hashing.h"

using namespace hashing;

hash_t hashing::in_range(hash_t v, hash_t lo, hash_t hi) {
    REQUIRE_GT(hi, lo, "upper hash value does not exceed lower value");
    v = fnv_hash(v + 312194ul);
    v %= hi - lo;
    return v + lo;
}

hash_t hashing::in_range(const std::vector<hash_t> &v, hash_t lo, hash_t hi) {
    REQUIRE_FALSE(v.empty(), "there must be as least one hashing value");
    auto out = in_range(v[0], lo, hi);
    for (size_t i = 1ul; i < v.size(); ++i) out = in_range(out + 4321 * v[i], lo, hi);
    return out;
}

std::vector<hash_t> hashing::in_range(const std::vector<hash_t> &v, size_t ngen, hash_t lo, hash_t hi, bool sorted) {
    std::vector<hash_t> out;
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

std::vector<hash_t> hashing::in_range(hash_t v, size_t ngen, hash_t lo, hash_t hi, bool sorted) {
    return in_range(std::vector<hash_t>{v}, ngen, lo, hi, sorted);
}


std::vector<hash_t> hashing::unique_in_range(const std::vector<hash_t> &v, size_t ngen, hash_t lo, hash_t hi, bool sorted) {
    REQUIRE_LE(lo + ngen, hi, "number of unique values can't exceed the range of allowed values");
    std::vector<hash_t> out;
    out.reserve(ngen);
    std::set<hash_t> set;
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

std::vector<hash_t> hashing::unique_in_range(hash_t v, size_t ngen, hash_t lo, hash_t hi, bool sorted) {
    return unique_in_range(std::vector<hash_t>{v}, ngen, lo, hi, sorted);
}
