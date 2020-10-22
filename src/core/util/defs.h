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
#include <cstdlib>
#include "consts.h"

#ifndef NDEBUG
#ifdef VERBOSE
#define VERBOSE_DEBUGGING
#endif
#endif

#ifdef NDEBUG
#define	ASSERT(e) {}
#else
#define ASSERT(e) \
{ \
    if(!(e)){ \
        fputs("Assertion failed: \n\t", stderr); \
        fputs(#e, stderr); \
        fprintf(stderr, "\nin file \"%s\", line %d\n\n", __FILE__, __LINE__); \
        fflush(stderr); \
        abort(); \
    } \
}
#endif

#define HERE() std::cout << "file\""<<__FILE__<<"\", line " << __LINE__ << std::endl;

#ifdef NDEBUG
#define DBVAR(v)
#else
#define DBVAR(v) std::cout << "variable " << #v << " = " << v << " in " << "file\""<<__FILE__<<"\", line " << __LINE__ << std::endl;
#endif


namespace defs {
    const std::string assets_root = PROJECT_ROOT"/assets";
    typedef std::vector<size_t> inds;
    typedef int mpi_count;
    typedef std::vector<mpi_count> mpi_counts;
    /*
     * determinant decoding / sampling requires statically allocated arrays.
     * since the exact size of the determinant is not known at compile time,
     * this constant should be set to ensure plenty of space for the storage
     * of set/clear bit positions
     */
    constexpr size_t det_work_size = 512;
    typedef std::array<size_t, det_work_size> det_work;
    typedef std::pair<size_t, size_t> pair;
    //typedef std::complex<double> ham_t;
    typedef double ham_t;
    typedef ham_t wf_t;
    typedef typename consts::component_t<ham_t>::type ham_comp_t;
    typedef typename consts::component_t<wf_t>::type wf_comp_t;
    typedef double prob_t;
    typedef uint64_t hash_t;
    typedef uint64_t data_t;
    constexpr size_t nbyte_data = sizeof(data_t);
    constexpr size_t nbit_data = CHAR_BIT * nbyte_data;
    const size_t isym_1e = 2;
    const size_t isym_2e = 4;

    // width of the cache line in bytes
    constexpr size_t ndata_cacheline = 8;
    constexpr size_t ncacheline_byte = ndata_cacheline * sizeof(data_t);
}
#endif //M7_DEFS_H
