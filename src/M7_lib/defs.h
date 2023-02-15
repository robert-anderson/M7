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

/**
 * the integral type used throughout the code to count and number items. thus, it is unsigned
 */
typedef std::size_t uint_t;
/**
 * the vector STL container is used ubiquitously, so here a shortcut is created
 */
template<typename T>
using v_t = std::vector<T>;
/**
 * vector of such integers
 */
typedef v_t<uint_t> uintv_t;
/**
 * array of such integers
 */
template<uint_t n>
using uinta_t = std::array<uint_t, n>;
/**
 * pair of such integers
 */
typedef std::pair<uint_t, uint_t> uintp_t;

typedef std::string str_t;
typedef v_t<str_t> strv_t;
typedef std::pair<str_t, str_t> strp_t;
typedef std::map<str_t, str_t> strm_t;


#ifndef NDEBUG
constexpr bool c_enable_debug = true;
#else
constexpr bool c_enable_debug = false;
#endif

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

typedef double prob_t;
typedef unsigned char buf_t;
typedef unsigned char mae_ind_t;
typedef buf_t bos_occ_t;

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

constexpr uint_t c_mbf_type_ind = MBF_TYPE;

#ifdef ENABLE_BOSONS
constexpr bool c_enable_bosons = true;
#else
constexpr bool c_enable_bosons = false;
#endif

#ifdef ENABLE_FERMIONS
constexpr bool c_enable_fermions = true;
#else
constexpr bool c_enable_fermions = false;
#endif

//  nroot, nreplica
constexpr uint_t c_ndim_wf = 2;
//  nroot
constexpr uint_t c_ndim_root = 1;
typedef uinta_t<c_ndim_wf> wf_inds_t;

#ifdef ENABLE_TCHINT
constexpr bool c_enable_tchint = true;
#else
constexpr bool c_enable_tchint = false;
#endif

#ifdef ENABLE_LOCAL_LOGGING
constexpr bool c_enable_local_logging = true;
#else
constexpr bool c_enable_local_logging = false;
#endif

#endif //M7_DEFS_H
