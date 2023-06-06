//
// Created by rja on 24/03/23.
//

#include "Stats.h"

wf::Stats::Stats(const NdFormat<c_ndim_wf>& format, uint_t nshift_space) :
        m_ninitiator(format), m_nocc_mbf_by_shift_space({nshift_space}),
        m_nw(format), m_nw_by_shift_space(format + nshift_space), m_l2_norm_square(format),
        m_nspawned(format), m_nannihilated(format) {
    m_summed = {&m_ninitiator, &m_nocc_mbf_by_shift_space, &m_nocc_mbf, &m_nw, &m_nw_by_shift_space,
                &m_l2_norm_square, &m_nspawned, &m_nannihilated, &m_nlarge_ci};
}