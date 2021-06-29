//
// Created by rja on 18/05/2020.
//

#ifndef M7_MAGNITUDELOGGER_H
#define M7_MAGNITUDELOGGER_H


#include <cstddef>
#include <src/defs.h>
#include <src/core/parallel/Reducible.h>
#include <src/core/parallel/Epoch.h>
#include "src/core/config/FciqmcConfig.h"

class MagnitudeLogger {
    const fciqmc_config::Propagator &m_opts;
    size_t m_nsingle = 0;
    size_t m_ndouble = 0;
    // highest magnitudes
    Reducible<defs::ham_comp_t> m_hi_mag_single;
    Reducible<defs::ham_comp_t> m_hi_mag_double;

    Epoch m_enough_singles_for_dynamic_tau;
    Epoch m_enough_doubles_for_dynamic_tau;

    defs::prob_t psingle_guess(size_t nsite, size_t nelec){
        return 1.0/(integer_utils::combinatorial(2*nsite-nelec, 2)/double(2*nsite-nelec)+1);
    }

public:
    // the recommended timestep based on the hi_mag and the maximum acceptable bloom
    defs::prob_t m_psingle;
    double m_tau;

    MagnitudeLogger(const fciqmc_config::Propagator &opts, defs::prob_t m_psingle);

    MagnitudeLogger(const fciqmc_config::Propagator &opts, size_t nsite, size_t nelec);

    void log(size_t nexcit, defs::ham_t helem, defs::prob_t prob);

    void synchronize(size_t icycle);

};


#endif //M7_MAGNITUDELOGGER_H
