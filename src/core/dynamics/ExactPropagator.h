//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include "src/core/config/FciqmcConfig.h"
#include "src/core/excititer/ExcitIterGroup.h"
#include "Propagator.h"

class ExactPropagator : public Propagator {
    /**
     * sends the generated excitations even if the corresponding hamiltonian matrix element is zero. Useful for testing
     * rank-2 RDMs since these spawns will make the exact contributions.
     */
    const bool m_only_nonzero_h_spawns;
    ExcitIterGroup m_excit_iters;
    MagnitudeLogger m_mag_log;

public:
    ExactPropagator(const Hamiltonian &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt,
                    bool only_nonzero_h_spawns=true);

    void diagonal(Wavefunction &wf, const size_t& ipart) override;

    void off_diagonal(Wavefunction &wf, const size_t& ipart) override;

    void update(const size_t &icycle, const Wavefunction &wf) override;
};

#endif //M7_EXACTPROPAGATOR_H
