//
// Created by rja on 13/04/2021.
//

#ifndef M7_FCIQMCSTATS_H
#define M7_FCIQMCSTATS_H

#include <src/core/dynamics/Propagator.h>
#include "src/core/field/Row.h"
#include "StatsTable.h"

struct FciqmcStatsRow : Row {

    NdFormat<defs::ndim_wf> m_wf_format;

    field::Number<size_t> m_icycle;
    field::Number<double> m_tau;
    field::Numbers<defs::ham_comp_t, defs::ndim_wf> m_shift;
    field::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    field::Numbers<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    field::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_spawned;
    field::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_annihilated;
    field::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy_num;
    field::Numbers<defs::wf_t, defs::ndim_wf> m_ref_weight;
    field::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy;
    field::Numbers<defs::ham_comp_t, defs::ndim_wf> m_l2_norm;
    field::Numbers<size_t, defs::ndim_wf> m_ninitiator;
    field::Numbers<size_t, defs::ndim_wf> m_nocc_mbf;
    field::Numbers<int, defs::ndim_wf> m_delta_nocc_mbf;
    field::Numbers<defs::prob_t, 1ul> m_exlvl_probs;
    field::Numbers<defs::ham_t, defs::ndim_wf> m_uniform_twf_num;
    field::Numbers<defs::ham_t, defs::ndim_wf> m_weighted_twf_num;
    field::Numbers<defs::ham_t, defs::ndim_wf> m_weighted_twf_denom;
    field::Numbers<defs::ham_comp_t, defs::ndim_wf> m_reweighting_factor;

    FciqmcStatsRow(Propagator& prop);
};

typedef StatsTable<FciqmcStatsRow> FciqmcStats;

#endif //M7_FCIQMCSTATS_H
