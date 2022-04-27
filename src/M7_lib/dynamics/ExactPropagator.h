//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include <M7_lib/config/FciqmcConfig.h>

#include "Propagator.h"
#include "M7_lib/foreach/ConnForeachGroup.h"

class ExactPropagator : public Propagator {
    /**
     * sends the generated excitations even if the corresponding hamiltonian matrix element is zero. Useful for testing
     * rank-2 RDMs since these spawns will make the exact contributions.
     */
    const bool m_only_nonzero_h_spawns;
    ConnForeachGroup m_conn_iters;
    MagnitudeLogger m_mag_log;

public:
    ExactPropagator(const Hamiltonian &ham, const conf::Document &opts, const Wavefunction& wf,
                    bool only_nonzero_h_spawns=true);

    void diagonal(Wavefunction &wf, const size_t& ipart) override;

    void off_diagonal(Wavefunction &wf, const size_t& dst_mbf) override;

    void update(const size_t &icycle, const Wavefunction &wf) override;
};

#endif //M7_EXACTPROPAGATOR_H
