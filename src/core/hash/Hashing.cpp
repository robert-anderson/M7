//
// Created by anderson on 1/7/22.
//

#include "Hashing.h"

defs::hash_t hashing::in_range(defs::hash_t v, const defs::hash_t &lo, const defs::hash_t &hi) {
    ASSERT(hi>lo);
    v = fnv_hash(v+312194ul);
    v%=hi-lo;
    return v+lo;
}

defs::hash_t hashing::in_range(const std::vector<defs::hash_t> &v, const defs::hash_t &lo, const defs::hash_t &hi) {
    ASSERT(!v.empty());
    auto out = in_range(v[0], lo, hi);
    for (size_t i=1ul; i<v.size(); ++i) out = in_range(out+4321*v[i], lo, hi);
    return out;
}

std::vector<defs::hash_t>
hashing::in_range(const std::vector<defs::hash_t> &v, size_t ngen, const defs::hash_t &lo, const defs::hash_t &hi,
                  bool sorted) {
    ASSERT(lo+ngen<=hi);
    std::vector<defs::hash_t> out;
    out.reserve(ngen);
    std::set<defs::hash_t> set;
    auto vatt = v;
    vatt.push_back(0);
    while (set.size()!=ngen){
        auto r = in_range(vatt, lo, hi);
        while (set.find(r)!=set.end()) if (++r==hi) r = lo;
        out.push_back(r);
        set.insert(r);
        ++vatt.back();
    }
    if (sorted) std::sort(out.begin(), out.end());
    return out;
}

std::vector<defs::hash_t>
hashing::in_range(defs::hash_t v, size_t ngen, const defs::hash_t &lo, const defs::hash_t &hi, bool sorted) {
    return in_range(std::vector<defs::hash_t>{v}, ngen, lo, hi, sorted);
}
