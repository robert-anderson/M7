//
// Created by Robert John Anderson on 2020-01-05.
//

#ifndef M7_DEFS_H
#define M7_DEFS_H

#include <string>
#include <map>
#include <complex>
#include <iostream>
#include <vector>
#include <queue>
#include <exception>
#include <array>
#include <climits>
#include <cstdlib>
#include <M7_lib/util/consts.h>

/**
 * basic debugging macro for verification of low-level functionality, in general prefer the macros MpiAssert.h defines
 */
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
        std::abort(); \
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

#ifndef NDEBUG
    constexpr bool enable_debug = true;
#else
    constexpr bool enable_debug = false;
#endif
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
    typedef std::pair<size_t, size_t> pair;
#ifdef ENABLE_COMPLEX
    constexpr bool enable_complex = true;
#else
    constexpr bool enable_complex = false;
#endif

    typedef double ham_comp_t;
    typedef ham_comp_t wf_comp_t;
    typedef std::conditional<enable_complex, std::complex<ham_comp_t>, ham_comp_t>::type ham_t;
    typedef std::map<std::string, std::string> info_map_t;
    typedef ham_t wf_t;
    typedef double prob_t;
    typedef uint64_t hash_t;
    typedef char buf_t;
    typedef unsigned char mev_ind_t;
    typedef unsigned char bos_occ_t;
    constexpr int undefined_ms2 = std::numeric_limits<int>::max();
    constexpr size_t max_bos_occ = std::numeric_limits<bos_occ_t>::max();
    constexpr size_t nbyte_word = sizeof(size_t);
    constexpr size_t nbit_word = CHAR_BIT * nbyte_word;
    /**
     * Hamiltonian-parametrising coefficients with magnitude below this value are considered to be zero, and wherever
     * equality between such values is to be asserted, this same tolerance is applied
     */
    constexpr ham_comp_t helem_tol = 1e-9;

    /**
     * "exsigs", short for "excitation signatures", encode in a single word the number of each type of second-quantised
     * operator in an operator product.
     */
    /**
     * number of bits in the signature representing each number of fermion SQ operators
     */
    static constexpr size_t nbit_exsig_nop_frm = 3;
    /**
     * number of bits in the signature representing each number of boson SQ operators
     */
    static constexpr size_t nbit_exsig_nop_bos = 2;
    /**
     * mask and max value for extraction of a number of fermion SQ operators
     */
    static constexpr size_t exsig_nop_mask_frm = (1 << nbit_exsig_nop_frm)-1;
    /**
     * mask and max value for extraction of a number of boson SQ operators
     */
    static constexpr size_t exsig_nop_mask_bos = (1 << nbit_exsig_nop_bos)-1;
    /**
     * total number of distinct excitation signatures that can be stored
     */
    static constexpr size_t nexsig = (1 << (2 * nbit_exsig_nop_frm + 2*nbit_exsig_nop_bos));

    enum MbfTypeInd {Frm, FrmBos, Bos};
    /**
     * Many-body basis function definitions
     *  0: fermion only (determinant basis)
     *  1: fermion-boson (determinant-permanent product basis)
     *  2: boson only (permanent basis)
     *  3: fermion only spin-adapted (CSF basis) (TODO: not currently implemented)
     */
#ifndef MBF_TYPE
    #define MBF_TYPE 0
#endif

#if (MBF_TYPE==1) || (MBF_TYPE==2)
    #define ENABLE_BOSONS 1
#endif

#if (MBF_TYPE==0) || (MBF_TYPE==1)
    #define ENABLE_FERMIONS
#endif

    constexpr size_t mbf_type_ind = MBF_TYPE;

#ifdef ENABLE_BOSONS
    constexpr bool enable_bosons = true;
#else
    constexpr bool enable_bosons = false;
#endif

#ifdef ENABLE_FERMIONS
    constexpr bool enable_fermions = true;
#else
    constexpr bool enable_fermions = false;
#endif

    //  nroot, nreplica
    constexpr size_t ndim_wf = 2;
    //  nroot
    constexpr size_t ndim_root = 1;
    typedef std::array<size_t, ndim_wf> wf_inds_t;
}
#endif //M7_DEFS_H
