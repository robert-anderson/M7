//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_FCIQMCSTATS_H
#define M7_FCIQMCSTATS_H

#include <M7_lib/dynamics/Propagator.h>
#include <M7_lib/field/Row.h>

#include "StatsTable.h"

struct FciqmcStatsRow : StatsRow {

    NdFormat<defs::ndim_wf> m_wf_format;

    statistic::Number<size_t> m_icycle;
    statistic::Number<double> m_tau;
    statistic::Numbers<defs::ham_comp_t, defs::ndim_wf> m_shift;
    statistic::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    statistic::Numbers<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    statistic::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_spawned;
    statistic::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_annihilated;
    statistic::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy_num;
    statistic::Numbers<defs::wf_t, defs::ndim_wf> m_ref_weight;
    statistic::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy;
    statistic::Numbers<defs::ham_comp_t, defs::ndim_wf> m_l2_norm;
    statistic::Numbers<size_t, defs::ndim_wf> m_ninitiator;
    statistic::Numbers<size_t, defs::ndim_wf> m_nocc_mbf;
    statistic::Numbers<int, defs::ndim_wf> m_delta_nocc_mbf;
    statistic::Numbers<defs::prob_t, 1ul> m_exlvl_probs;
    statistic::Numbers<defs::ham_comp_t, defs::ndim_wf> m_reweighting_factor;

    FciqmcStatsRow(Propagator& prop);
};

typedef StatsTable<FciqmcStatsRow> FciqmcStats;

#endif //M7_FCIQMCSTATS_H
