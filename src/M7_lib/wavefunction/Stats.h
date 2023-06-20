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
         * number of MBFs with any associated weight in any part
         */
        reduction::cyclic::Scalar<int64_t, false> m_nocc_mbf;
        /**
         * number of MBFs with any associated weight for each shift space
         */
        reduction::cyclic::NdArray<int64_t, 1, false> m_nocc_mbf_by_shift_space;
        /**
         * L1 norm of each part of the WF
         */
        reduction::cyclic::NdArray<wf_comp_t, c_ndim_wf> m_nw;
        /**
         * L1 norm of each shift space, and each part of the WF
         */
        reduction::cyclic::NdArray<wf_comp_t, 1+c_ndim_wf> m_nw_by_shift_space;
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
        /**
         * number of MBFs in the large CI set
         */
        reduction::cyclic::Scalar<uint_t> m_nlarge_ci;

        Stats(const NdFormat<c_ndim_wf>& format, uint_t nshift_space);

        void all_sum() {
            reduction::all_sum(m_summed);
        }
    };
}


#endif //M7_WF_STATS_H
