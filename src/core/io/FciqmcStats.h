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
    fields::Numbers<defs::wf_comp_t, defs::ndim_wf> m_nwalker_annihilated;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_ref_proj_energy_num;
    fields::Numbers<defs::wf_t, defs::ndim_wf> m_ref_weight;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_ref_proj_energy;
    fields::Numbers<defs::ham_comp_t, defs::ndim_wf> m_l2_norm;
    fields::Numbers<size_t, defs::ndim_wf> m_ninitiator;
    fields::Numbers<size_t, defs::ndim_wf> m_nocc_onv;
    fields::Numbers<int, defs::ndim_wf> m_delta_nocc_onv;
    fields::Numbers<defs::prob_t, defs::ndim_wf> m_psingle;
    fields::Number<double> m_total_synchronization_overhead;
    fields::Number<double> m_propagate_loop_time;
    fields::Number<double> m_communication_time;
    fields::Number<double> m_annihilation_loop_time;
    fields::Number<double> m_total_cycle_time;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_uniform_twf_num;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_hubbard_twf_num;
    fields::Numbers<defs::ham_t, defs::ndim_wf> m_hubbard_twf_denom;

    FciqmcStatsRow(NdFormat<2> format) :
    m_wf_format(format),
    m_icycle(this, "Cycle number"),
    m_tau(this, m_wf_format, "Timestep"),
    m_shift(this, m_wf_format, "Diagonal shift"),
    m_nwalker(this, m_wf_format, "WF L1 norm (number of walkers)"),
    m_delta_nwalker(this, m_wf_format, "Walkers added this cycle"),
    m_nwalker_annihilated(this, m_wf_format, "Walkers annihilated this cycle"),
    m_ref_proj_energy_num(this, m_wf_format, "Numerator of reference-projected energy estimator"),
    m_ref_weight(this, m_wf_format, "Reference weight"),
    m_ref_proj_energy(this, m_wf_format, "Reference-projected energy"),
    m_l2_norm(this, m_wf_format, "L2 norm of the wavefunction"),
    m_ninitiator(this, m_wf_format, "Initiator ONVs"),
    m_nocc_onv(this, m_wf_format, "Occupied ONVs"),
    m_delta_nocc_onv(this, m_wf_format, "Change in number of occupied ONVs"),
    m_psingle(this, m_wf_format, "Probability of attempting to draw a single excitation"),
    m_total_synchronization_overhead(this, "Total time waited at MPI_Barrier"),
    m_propagate_loop_time(this, "Time spent in loop over occupied ONVs"),
    m_communication_time(this, "Time spent in communicating spawns"),
    m_annihilation_loop_time(this, "Time spent in annihilation loop"),
    m_total_cycle_time(this, "Total cycle time"),
    m_uniform_twf_num(this, m_wf_format, "Numerator of uniform TWF-projected energy estimator"),
    m_hubbard_twf_num(this, m_wf_format, "Numerator of Hubbard-Holstein-like TWF-projected energy estimator"),
    m_hubbard_twf_denom(this, m_wf_format, "Denominator of Hubbard-Holstein-like TWF-projected energy estimator")
    {}

    FciqmcStatsRow(size_t nroot, size_t nreplica) :
            FciqmcStatsRow({{nroot, nreplica}, {"nroot", "nreplica"}}){}
};

typedef StatsTable<FciqmcStatsRow> FciqmcStats;

#endif //M7_FCIQMCSTATS_H
