//
// Created by Robert John Anderson on 2020-01-05.
//

#ifndef M7_DEFS_H
#define M7_DEFS_H

#include <string>
#include <map>
#include <complex>
#include <vector>
#include <climits>
#include <limits>

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

namespace defs {

#ifndef NDEBUG
    constexpr bool c_enable_debug = true;
#else
    constexpr bool c_enable_debug = false;
#endif
    /**
     * the integral type used throughout the code to count and number items. thus, it is unsigned
     */
    typedef uint64_t i_t;
    /**
     * vector of such integers
     */
    typedef std::vector<i_t> ivec_t;
    /**
     * pair of such integers
     */
    typedef std::pair<i_t, i_t> ipair_t;

#ifdef ENABLE_COMPLEX_HAM
    constexpr bool c_enable_complex_ham = true;
#else
    constexpr bool c_enable_complex_ham = false;
#endif

#ifdef ENABLE_COMPLEX_WF
    constexpr bool c_enable_complex_wf = true;
#else
    constexpr bool c_enable_complex_wf = false;
#endif

    // type of each arithmetic component of the Hamiltonian matrix elements
    typedef double ham_comp_t;
    // overall type of the Hamiltonian matrix elements (taking possible complex enabling into account)
    typedef std::conditional<c_enable_complex_ham, std::complex<ham_comp_t>, ham_comp_t>::type ham_t;

    // type of each arithmetic component of the wavefunctions (also MAEs and spawns emitted by the propagators)
    typedef ham_comp_t wf_comp_t;
    // overall type of the wavefunction elements (taking possible complex enabling into account)
    typedef std::conditional<c_enable_complex_wf, std::complex<wf_comp_t>, wf_comp_t>::type wf_t;

    typedef std::map<std::string, std::string> info_map_t;
    typedef double prob_t;
    typedef char buf_t;
    typedef unsigned char mev_ind_t;
    typedef unsigned char bos_occ_t;

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
