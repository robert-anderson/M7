//
// Created by rja on 25/08/2021.
//

#ifndef M7_EXCITITERGROUP_H
#define M7_EXCITITERGROUP_H

#include "Hubbard1dSingles.h"

using namespace exsig_utils;
using namespace excititers;

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

    void add(std::unique_ptr<ExcitIter> &&excit_iter, size_t exsig);

public:
    /**
     * infer from the given Hamiltonian exactly which excitation iterators are required, and initialize them
     * @param ham
     *  general Hamiltonian object whose excitation level information is queried to determine the required exsig-specific
     *  excitation iterator objects
     */
    ExcitIterGroup(const Hamiltonian &ham) {
        if (ham.m_frm.m_model_attrs.is_hubbard_1d() || ham.m_frm.m_model_attrs.is_hubbard_1d_pbc()) {
            add(std::unique_ptr<ExcitIter>(new Hubbard1dSingles(ham)), ex_single);
        } else {
            add(std::unique_ptr<ExcitIter>(new FrmConserve(ham, ex_single)), ex_single);
        }
        if (ham.m_frm.m_contribs_2200.is_nonzero(ex_double)) {
            add(std::unique_ptr<ExcitIter>(new FrmConserve(ham, ex_double)), ex_double);
        }
        if (ham.m_bos.m_nboson_max) {
//        m_exgens[conn_utils::encode_exsig(0, 0, 1, 0)] =
//                std::unique_ptr<ExcitGen>(new UniformHolstein(ham, prng, true));
//        m_exgens[conn_utils::encode_exsig(0, 0, 0, 1)] =
//                std::unique_ptr<ExcitGen>(new UniformHolstein(ham, prng, false));
        }

        init();
    }
};


#endif //M7_EXCITITERGROUP_H
