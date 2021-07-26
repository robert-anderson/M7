//
// Created by rja on 23/07/2021.
//

#include "FciqmcStats.h"

FciqmcStatsRow::FciqmcStatsRow(NdFormat<2> format) :
        m_wf_format(format),
        m_icycle(this, "Cycle number"),
        m_tau(this, m_wf_format, "Timestep"),
        m_shift(this, m_wf_format, "Diagonal shift"),
        m_nwalker(this, m_wf_format, "WF L1 norm (number of walkers)"),
        m_delta_nwalker(this, m_wf_format, "Walkers added this cycle"),
        m_nwalker_spawned(this, m_wf_format, "Walkers spawned this cycle"),
        m_nwalker_annihilated(this, m_wf_format, "Walkers annihilated this cycle"),
        m_ref_proj_energy_num(this, m_wf_format, "Numerator of reference-projected energy estimator"),
        m_ref_weight(this, m_wf_format, "Reference weight"),
        m_ref_proj_energy(this, m_wf_format, "Reference-projected energy"),
        m_l2_norm(this, m_wf_format, "L2 norm of the wavefunction"),
        m_ninitiator(this, m_wf_format, "Initiator ONVs"),
        m_nocc_mbf(this, m_wf_format, "Occupied ONVs"),
        m_delta_nocc_mbf(this, m_wf_format, "Change in number of occupied ONVs"),
        m_psingle(this, m_wf_format, "Probability of attempting to draw a single excitation"),
        m_uniform_twf_num(this, m_wf_format, "Numerator of uniform TWF-projected energy estimator"),
        m_weighted_twf_num(this, m_wf_format, "Numerator of weighted TWF-projected energy estimator"),
        m_weighted_twf_denom(this, m_wf_format, "Denominator of weighted TWF-projected energy estimator"),
        m_reweighting_factor(this, m_wf_format, "Reweighting factor for population "
                                                "control unbiasing")
{}

FciqmcStatsRow::FciqmcStatsRow(size_t nroot, size_t nreplica) :
        FciqmcStatsRow({{nroot, nreplica}, {"nroot", "nreplica"}}){}
