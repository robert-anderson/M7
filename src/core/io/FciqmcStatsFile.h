//
// Created by rja on 28/02/2020.
//

#ifndef M7_FCIQMCSTATSFILE_H
#define M7_FCIQMCSTATSFILE_H

#include "StatsFile.h"
#include "Options.h"
struct FciqmcStatsFile : public StatsFile {
public:
    StatsField<size_t> m_cycle_number;
    StatsField<defs::ham_comp_t> m_diagonal_shift;
    StatsField<double> m_timestep;
    StatsField<defs::ham_t> m_ref_proj_energy_num;
    StatsField<defs::ham_t> m_ref_weight;
    StatsField<defs::ham_t> m_ref_proj_energy;
    StatsField<defs::ham_comp_t> m_nwalker;
    StatsField<defs::ham_comp_t> m_nw_growth_rate;
    StatsField<defs::wf_comp_t> m_nw_at_doubles;
    StatsField<defs::wf_comp_t> m_ref_candidate_weight;
    StatsField<size_t> m_ninitiator;
    StatsField<defs::wf_t> m_aborted_weight;
    StatsField<size_t> m_noccupied_det;
    StatsField<defs::prob_t> m_psingle;
    StatsField<double> m_iter_time;
    StatsField<double> m_prop_time;
    StatsField<double> m_comm_time;
    StatsField<double> m_anni_time;

    explicit FciqmcStatsFile(const Options &input);
};


#endif //M7_FCIQMCSTATSFILE_H
