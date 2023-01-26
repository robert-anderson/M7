//
// Created by Robert J. Anderson on 27/02/2020.
//

#ifndef M7_EXACT_LINEAR_PROPAGATOR_H
#define M7_EXACT_LINEAR_PROPAGATOR_H

#include <M7_lib/conf/Conf.h>

#include "M7_lib/propagator/Propagator.h"
#include "M7_lib/foreach/ConnForeachGroup.h"

class ExactLinear : public Propagator {
    /**
     * sends the generated excitations even if the corresponding hamiltonian matrix element is zero. Useful for testing
     * rank-2 RDMs since these spawns will make the exact contributions.
     */
    const bool m_only_nonzero_h_spawns;
    ConnForeachGroup m_conn_iters;
    MagnitudeLogger m_mag_log;

public:
    ExactLinear(const Hamiltonian& ham, const conf::Document& opts, const Wavefunction& wf,
                bool only_nonzero_h_spawns=true);

    void diagonal(Wavefunction &wf, Walker& walker, uint_t ipart) override;

    void off_diagonal(Wavefunction &wf, const Walker& walker, uint_t dst_mbf, bool initiator) override;

    void update(uint_t icycle, const Wavefunction &wf) override;
};

#endif //M7_EXACT_LINEAR_PROPAGATOR_H
