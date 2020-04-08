//
// Created by rja on 28/02/2020.
//

#ifndef M7_FCIQMCSTATSFILE_H
#define M7_FCIQMCSTATSFILE_H

#include "StatsFile.h"
#include "InputOptions.h"
struct FciqmcStatsFile : public StatsFile {
public:
    StatsField<size_t> m_cycle_number;
    StatsField<defs::ham_comp_t> m_diagonal_shift;
    StatsField<double> m_timestep;
    StatsField<defs::ham_t> m_ref_proj_energy_num;
    StatsField<defs::ham_t> m_ref_weight;
    StatsField<defs::ham_comp_t> m_nwalker;
    StatsField<size_t> m_ninitiator;
    StatsField<defs::ham_comp_t> m_aborted_weight;
    StatsField<size_t> m_noccupied_det;

    explicit FciqmcStatsFile(const InputOptions &input):
    StatsFile(input.stats_path),
    m_cycle_number(this, 1, "Cycle number"),
    m_diagonal_shift(this, 1, "Diagonal shift"),
    m_timestep(this, 1, "Timestep"),
    m_ref_proj_energy_num(this, 1, "Reference projected energy numerator"),
    m_ref_weight(this, 1, "Reference weight"),
    m_nwalker(this, 1, "Total number of walkers"),
    m_ninitiator(this, 1, "Number of initiator determinants"),
    m_aborted_weight(this, 1, "Aborted weight"),
    m_noccupied_det(this, 1, "Number of occupied determinants"){}
};


#endif //M7_FCIQMCSTATSFILE_H
