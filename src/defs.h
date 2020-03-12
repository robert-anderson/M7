//
// Created by Robert John Anderson on 2020-01-05.
//

#ifndef M7_DEFS_H
#define M7_DEFS_H

#include <string>
#include <complex>
#include <vector>
#include "consts.h"

namespace defs {
    const std::string assets_root = "../assets";
    typedef std::vector<size_t> inds;
    typedef std::pair<size_t, size_t> pair;
    typedef std::complex<double> ham_t;
    typedef typename consts::component_t<ham_t>::type ham_comp_t;
    typedef double prob_t;
    typedef uint64_t data_t;
    const size_t isym_1e = 2;
    const size_t isym_2e = 4;

    // width of the cache line in bytes
    const size_t cache_line_size = 64;
}

#endif //M7_DEFS_H
