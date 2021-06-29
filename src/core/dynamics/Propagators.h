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

    static std::unique_ptr<Propagator> get(const Hamiltonian<>& ham, const fciqmc_config::Document& opts,
                                           const NdFormat<defs::ndim_wf>& wf_fmt){
        if (opts.m_propagator.m_exact)
            return std::unique_ptr<Exact>(new Exact{ham, opts, wf_fmt});
        else
            return std::unique_ptr<Stoch>(new Stoch{ham, opts, wf_fmt});
    }

}

#endif //M7_PROPAGATORS_H
