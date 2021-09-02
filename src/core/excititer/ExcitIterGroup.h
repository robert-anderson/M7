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
    ExcitIterGroup(const Hamiltonian &ham) {

        bool any_singles =
                ham.m_frm.m_contribs_1100.is_nonzero(ex_single) || ham.m_frm.m_contribs_2200.is_nonzero(ex_single);
        if (any_singles) {
            if (ham.m_frm.m_model_attrs.is_hubbard_1d() || ham.m_frm.m_model_attrs.is_hubbard_1d_pbc()) {
                add(std::unique_ptr<ExcitIter>(new excititers::Hubbard1dSingles(ham)));
            } else {
                add(std::unique_ptr<ExcitIter>(new excititers::FrmConserve(ham, ex_single)));
            }
        }
        if (ham.m_frm.m_contribs_2200.is_nonzero(ex_double)) {
            add(std::unique_ptr<ExcitIter>(new excititers::FrmConserve(ham, ex_double)));
        }

        if (ham.m_bos.m_nboson_max) {
            bool any_pures;
            any_pures = ham.m_ladder.m_contribs_0010.is_nonzero(exsig_utils::ex_0010)
                             || ham.m_ladder.m_contribs_1110.is_nonzero(exsig_utils::ex_0010);
            if (any_pures) {
                if (ham.m_ladder.is_holstein())
                    add(std::unique_ptr<ExcitIter>(new excititers::LadderPureHolstein(ham, exsig_utils::ex_0010)));
                else
                    add(std::unique_ptr<ExcitIter>(new excititers::LadderPure(ham, exsig_utils::ex_0010)));
            }
            any_pures = ham.m_ladder.m_contribs_0001.is_nonzero(exsig_utils::ex_0001)
                             || ham.m_ladder.m_contribs_1101.is_nonzero(exsig_utils::ex_0001);
            if (any_pures) {
                if (ham.m_ladder.is_holstein())
                    add(std::unique_ptr<ExcitIter>(new excititers::LadderPureHolstein(ham, exsig_utils::ex_0001)));
                else
                    add(std::unique_ptr<ExcitIter>(new excititers::LadderPure(ham, exsig_utils::ex_0001)));
            }

            bool any_hopping;
            any_hopping = ham.m_ladder.m_contribs_1110.is_nonzero(exsig_utils::ex_1110);
            if (any_hopping)
                add(std::unique_ptr<ExcitIter>(new excititers::LadderHopping(ham, exsig_utils::ex_1110)));
            any_hopping = ham.m_ladder.m_contribs_1101.is_nonzero(exsig_utils::ex_1101);
            if (any_hopping)
                add(std::unique_ptr<ExcitIter>(new excititers::LadderHopping(ham, exsig_utils::ex_1101)));
        }

        init();
    }

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
