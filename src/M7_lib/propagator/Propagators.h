//
// Created by Robert J. Anderson on 10/06/2021.
//

#ifndef M7_PROPAGATORS_H
#define M7_PROPAGATORS_H

#include <memory>

#include "ExactLinear.h"
#include "StochLinear.h"
#include "M7_lib/util/Pointer.h"

namespace props {

    typedef ExactLinear Exact;
    typedef StochLinear Stoch;

    static std::unique_ptr<Propagator> get(const Hamiltonian &ham, const conf::Document &opts,
                                           const Wavefunction &wf) {
        /*
         * if the RDM contributions due to connections of the reference are not to be explicitly included on average,
         */
        bool only_nonzero_h_spawns = !opts.m_av_ests.m_rdm.m_enabled;
        if (!opts.m_propagator.m_stochastic)
            return ptr::smart::make_unique<Exact>(ham, opts, wf, only_nonzero_h_spawns);
        else
            return ptr::smart::make_unique<Stoch>(ham, opts, wf);
    }

}

#endif //M7_PROPAGATORS_H
