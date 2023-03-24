//
// Created by rja on 24/03/23.
//

#ifndef M7_WF_STATS_H
#define M7_WF_STATS_H

#include "M7_lib/parallel/Reduction.h"

namespace wf {
    /**
     * Stats about changes to the walker populations are kept separate to cut down on clutter in the Vectors class
     */
    struct Stats {
        /**
         * collection of all reductions which are summed at the end of every cycle
         */
        v_t<reduction::Base*> m_summed;

        /**
         * number of initiator MBFs in each part of the WF
         */
        reduction::NdArray<uint_t, c_ndim_wf> m_ninitiator;
        /**
         * number of initiator MBFs in each part of the WF due to permanitiator status
         */
        reduction::NdArray<uint_t, c_ndim_wf> m_ninitiator_perma;
        /**
         * number of MBFs with any associated weight in any part
         */
        reduction::cyclic::Scalar<int64_t, false> m_nocc_mbf;
        /**
         * L1 norm of each part of the WF
         */
        reduction::cyclic::NdArray<wf_comp_t, c_ndim_wf> m_nwalker;
        /**
         * square of the L2 norm of each part of the WF
         */
        reduction::cyclic::NdArray<wf_comp_t, c_ndim_wf> m_l2_norm_square;
        /**
         * number of walkers received in spawning process
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_nspawned;
        /**
         * number of walkers annihilated in the loop_over_spawned method for each part
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_nannihilated;

        Stats(const NdFormat<c_ndim_wf>& format);
    };
}


#endif //M7_WF_STATS_H
