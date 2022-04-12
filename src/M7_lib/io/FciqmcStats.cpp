//
// Created by rja on 23/07/2021.
//

#include "FciqmcStats.h"

FciqmcStatsRow::FciqmcStatsRow(Propagator& prop) :
        m_wf_format(prop.m_wf_fmt),
        m_icycle(this, "Cycle number", false),
        m_tau(this, "Timestep"),
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
        m_exlvl_probs(prop.ncase_excit_gen() ? this : nullptr,
                      {{prop.ncase_excit_gen()}, {"excitation generator index"}},
                      "Probability of attempting to draw excitation level"),
        m_reweighting_factor(this, m_wf_format, "Reweighting factor for population control unbiasing")
{}
