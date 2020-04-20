//
// Created by Robert John Anderson on 2020-01-05.
//

#ifndef M7_DEFS_H
#define M7_DEFS_H

#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <exception>
#include <array>
#include <climits>
#include "consts.h"

#ifdef NDEBUG
#define	assert(e)
#else
#define ASSERT(e) \
{if(!(e)){throw std::runtime_error(std::string("\nAssertion \"" #e "\" failed in file " __FILE__ )+" line: "+std::to_string( __LINE__ ));}}
#endif

namespace defs {
    const std::string assets_root = "../assets";
    typedef std::vector<size_t> inds;
    /*
     * determinant decoding / sampling requires statically allocated arrays.
     * since the exact size of the determinant is not known at compile time,
     * this constant should be set to ensure plenty of space for the storage
     * of set/clear bit positions
     */
    constexpr size_t det_work_size = 512;
    typedef std::array<size_t, det_work_size> det_work;
    typedef std::pair<size_t, size_t> pair;
    typedef std::complex<double> ham_t;
    //typedef double ham_t;
    typedef ham_t wf_t;
    typedef typename consts::component_t<ham_t>::type ham_comp_t;
    typedef typename consts::component_t<wf_t>::type wf_comp_t;
    typedef double prob_t;
    typedef uint64_t hash_t;
    typedef uint64_t data_t;
    constexpr size_t nbit_data = CHAR_BIT * sizeof(data_t);
    const size_t isym_1e = 2;
    const size_t isym_2e = 4;

    // width of the cache line in bytes
    constexpr size_t cache_line_size = 64;
}
#endif //M7_DEFS_H
