//
// Created by rja on 28/02/2020.
//

#ifndef M7_FCIQMCSTATSFILE_H
#define M7_FCIQMCSTATSFILE_H

#include "StatsFile.h"

struct FciqmcStatsSpecifier : StatsSpecifier {
    StatsColumn<size_t> m_icycle;
    StatsColumn<defs::ham_comp_t> m_tau;
    StatsColumn<defs::ham_comp_t> m_shift;
    StatsColumn<defs::wf_t> m_nwalker;
    StatsColumn<defs::wf_t> m_delta_nwalker;
    StatsColumn<defs::ham_t> m_ref_proj_energy_num;
    StatsColumn<defs::wf_t> m_ref_weight;
    StatsColumn<defs::ham_comp_t> m_ref_proj_energy;
    StatsColumn<defs::ham_comp_t> m_l2_norm;
    StatsColumn<size_t> m_ninitiator;
    StatsColumn<size_t> m_nocc_onv;
    StatsColumn<defs::prob_t> m_psingle;
    FciqmcStatsSpecifier() :
    StatsSpecifier("FCIQMC"),
    m_icycle(this, "Cycle number"),
    m_tau(this, "Timestep"),
    m_shift(this, "Diagonal shift"),
    m_nwalker(this, "Total number of walkers"),
    m_delta_nwalker(this, "Walkers added this cycle"),
    m_ref_proj_energy_num(this, "Numerator of reference-projected energy estimator"),
    m_ref_weight(this, "Reference weight"),
    m_ref_proj_energy(this, "Reference-projected energy"),
    m_l2_norm(this, "L2 norm of the wavefunction"),
    m_ninitiator(this, "Initiator ONVs"),
    m_nocc_onv(this, "Occupied ONVs"),
    m_psingle(this, "Probability of attempting to draw a single excitation")
    {}
};

#endif //M7_FCIQMCSTATSFILE_H
