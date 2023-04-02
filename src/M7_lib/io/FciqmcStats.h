//
// Created by Robert J. Anderson on 13/04/2021.
//

#ifndef M7_FCIQMCSTATS_H
#define M7_FCIQMCSTATS_H

#include <M7_lib/propagator/Propagator.h>
#include <M7_lib/field/Row.h>

#include "StatsTable.h"
#include "M7_lib/observables/InstEsts.h"

struct FciqmcStatsRow : StatsRow {

    NdFormat<c_ndim_wf> m_wf_format;
    // format for excitation-level resolved stats for each wavefunction part
    NdFormat<c_ndim_wf+1> m_wf_exlvl_format;

    statistic::Number<uint_t> m_icycle;
    statistic::Number<double> m_tau;
    statistic::Numbers<ham_comp_t, c_ndim_wf> m_shift;
    statistic::Numbers<wf_comp_t, c_ndim_wf> m_nwalker;
    statistic::Numbers<wf_comp_t, c_ndim_wf> m_delta_nwalker;
    statistic::Numbers<wf_comp_t, c_ndim_wf> m_nwalker_spawned;
    statistic::Numbers<wf_comp_t, c_ndim_wf> m_nwalker_annihilated;
    statistic::Numbers<ham_t, c_ndim_wf> m_ref_proj_energy_num;
    statistic::Numbers<wf_t, c_ndim_wf> m_ref_weight;
    statistic::Numbers<ham_t, c_ndim_wf> m_ref_proj_energy;
    statistic::Numbers<ham_comp_t, c_ndim_wf> m_l2_norm;
    statistic::Numbers<uint_t, c_ndim_wf> m_ninitiator;
    statistic::Numbers<uint_t, c_ndim_wf+1> m_ninitiator_by_exlvl;
    statistic::Numbers<uint_t, c_ndim_wf> m_nocc_pmntr;
    statistic::Numbers<uint_t, c_ndim_wf+1> m_nocc_pmntr_by_exlvl;
    statistic::Numbers<uint_t, c_ndim_wf> m_nocc_mbf;
    statistic::Numbers<int, c_ndim_wf> m_delta_nocc_mbf;
    statistic::Numbers<wf_t, c_ndim_wf> m_spin_square_num;
    statistic::Numbers<prob_t, 1ul> m_exlvl_probs;

    FciqmcStatsRow(const Propagator& prop, const InstEsts& inst_ests, bool exlvl_resolved);
};

typedef StatsTable<FciqmcStatsRow> FciqmcStats;

#endif //M7_FCIQMCSTATS_H
