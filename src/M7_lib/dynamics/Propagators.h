//
// Created by rja on 10/06/2021.
//

#ifndef M7_PROPAGATORS_H
#define M7_PROPAGATORS_H

#include <memory>

#include "ExactPropagator.h"
#include "StochasticPropagator.h"

namespace props {

    typedef ExactPropagator Exact;
    typedef StochasticPropagator Stoch;

    static std::unique_ptr<Propagator> get(const Hamiltonian &ham, const conf::Document &opts,
                                           const Wavefunction &wf) {
        /*
         * if the RDM contributions due to connections of the reference are not to be explicitly included on average,
         *
         */
        bool only_nonzero_h_spawns = opts.m_av_ests.m_rdm.m_explicit_ref_conns;
        if (!opts.m_propagator.m_stochastic)
            return std::unique_ptr<Exact>(new Exact{ham, opts, wf, only_nonzero_h_spawns});
        else
            return std::unique_ptr<Stoch>(new Stoch{ham, opts, wf});
    }

}

#endif //M7_PROPAGATORS_H
