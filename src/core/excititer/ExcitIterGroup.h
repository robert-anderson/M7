//
// Created by rja on 25/08/2021.
//

#ifndef M7_EXCITITERGROUP_H
#define M7_EXCITITERGROUP_H

#include "ExcitIters.h"
#include "BodyFnTypes.h"

using namespace exsig_utils;
using namespace body_fn_types;

struct ExcitIterGroup {
    std::array<std::unique_ptr<ExcitIter>, defs::nexsig> m_excit_iters;
    /**
     * vector storing all active exsigs consecutively (i.e. those elements of m_exiters which are non-null pointers to
     * ExcitIter objects)
     */
    defs::inds m_active_exsigs;
    /**
     * indices which point to positions in m_active_exsigs corresponding to purely fermionic excitations
     */
    defs::inds m_frm_inds;

    /**
     * initialize vectors of exsigs with allocated ExcitIter instances. Also initialize the vector m_frm_inds, which
     * identifies positions of purely fermionic excitation generators from the general vector
     */
    void init();

    void add(std::unique_ptr<ExcitIter> &&excit_iter);

public:
    /**
     * infer from the given Hamiltonian exactly which excitation iterators are required, and initialize them
     * @param ham
     *  general Hamiltonian object whose excitation level information is queried to determine the required exsig-specific
     *  excitation iterator objects
     */
    explicit ExcitIterGroup(const Hamiltonian &ham);

    void log_breakdown() const;


    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_c_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_cd_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_ch_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_cdh_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_d_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_h_t &body_fn) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn);
    }

    template<typename mbf_t>
    void foreach(const mbf_t &mbf, const fn_dh_t<mbf_t> &body_fn, bool nonzero_h_only) {
        for (const auto &exsig: m_active_exsigs) m_excit_iters[exsig]->foreach<mbf_t>(mbf, body_fn, nonzero_h_only);
    }
};


#endif //M7_EXCITITERGROUP_H