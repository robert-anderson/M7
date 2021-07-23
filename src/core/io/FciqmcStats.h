//
// Created by rja on 13/04/2021.
//

#ifndef M7_FCIQMCSTATS_H
#define M7_FCIQMCSTATS_H

#include "src/core/field/Row.h"
#include "StatsTable.h"

struct FciqmcStatsRow : Row {

    NdFormat<defs::ndim_wf> m_wf_format;

    fields::Number<size_t> m_icycle;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_tau;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_shift;
    fields::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    fields::Numbers<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    fields::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_spawned;
    fields::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_annihilated;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy_num;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_ref_weight;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_l2_norm;
    fields::Numbers<size_t, defs::ndim_wf> m_ninitiator;
    fields::Numbers<size_t, defs::ndim_wf> m_nocc_onv;
    fields::Numbers<int, defs::ndim_wf> m_delta_nocc_onv;
    fields::Numbers<defs::prob_t, defs::ndim_wf> m_psingle;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_uniform_twf_num;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_weighted_twf_num;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_weighted_twf_denom;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_reweighting_factor;

    FciqmcStatsRow(NdFormat<2> format);

    FciqmcStatsRow(size_t nroot, size_t nreplica);
};

typedef StatsTable<FciqmcStatsRow> FciqmcStats;

#endif //M7_FCIQMCSTATS_H
