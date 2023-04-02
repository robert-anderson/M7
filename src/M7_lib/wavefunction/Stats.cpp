//
// Created by rja on 24/03/23.
//

#include "Stats.h"

wf::Stats::Stats(const NdFormat<c_ndim_wf>& format) :
        m_ninitiator(format), m_nocc_pmntr(format), m_nwalker(format),
        m_l2_norm_square(format), m_nspawned(format), m_nannihilated(format) {
    m_summed = {&m_ninitiator, &m_nocc_pmntr, &m_nocc_mbf, &m_nwalker, &m_l2_norm_square, &m_nspawned, &m_nannihilated};
}