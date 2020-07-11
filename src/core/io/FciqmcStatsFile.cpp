//
// Created by rja on 28/02/2020.
//

#include "FciqmcStatsFile.h"

FciqmcStatsFile::FciqmcStatsFile(const Options &input) :
        StatsFile(input.stats_path),
        m_cycle_number(this, 1, "Cycle number"),
        m_diagonal_shift(this, 1, "Diagonal shift", 12),
        m_timestep(this, 1, "Timestep"),
        m_ref_proj_energy_num(this, 1, "Reference projected energy numerator", 12, true),
        m_ref_weight(this, 1, "Reference weight", 12, true),
        m_ref_proj_energy(this, 1, "Reference projected energy", 12),
        m_nwalker(this, 1, "Total number of walkers"),
        m_nw_growth_rate(this, 1, "Walker growth rate"),
        m_ninitiator(this, 1, "Number of initiator determinants"),
        m_aborted_weight(this, 1, "Aborted weight"),
        m_noccupied_det(this, 1, "Number of occupied determinants"),
        m_psingle(this, 1, "Probability that a stochastic single will be attempted"),
        m_iter_time(this, 1, "Total iteration time (seconds)"),
        m_prop_time(this, 1, "Propagation loop time (seconds)"),
        m_comm_time(this, 1, "Walker communication time (seconds)"),
        m_anni_time(this, 1, "Annihilation time (seconds)")
        {
    if (!mpi::i_am_root())
        throw std::runtime_error(
                "FciqmcStatsFiles must only be instantiated on the root process");
}
