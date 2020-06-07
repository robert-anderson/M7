//
// Created by rja on 18/05/2020.
//

#ifndef M7_MAGNITUDELOGGER_H
#define M7_MAGNITUDELOGGER_H


#include <cstddef>
#include <src/core/util/defs.h>
#include <src/core/thread/PrivateStore.h>
#include <src/core/thread/Reduction.h>
#include <src/core/parallel/MPIWrapper.h>
#include <src/core/io/Options.h>

class MagnitudeLogger {
    const Options &m_input;
    size_t m_nsingle = 0;
    size_t m_ndouble = 0;
    PrivateStore<size_t> m_priv_nsingle;
    PrivateStore<size_t> m_priv_ndouble;
    // highest magnitudes
    defs::ham_comp_t m_hi_mag_single = 0;
    defs::ham_comp_t m_hi_mag_double = 0;

    PrivateStore<defs::ham_comp_t> m_priv_hi_mag_single;
    PrivateStore<defs::ham_comp_t> m_priv_hi_mag_double;

    bool m_enough_singles_for_dynamic_tau = false;
    bool m_enough_doubles_for_dynamic_tau = false;

public:
    // the recommended timestep based on the hi_mag and the maximum acceptable bloom
    double m_tau;
    defs::prob_t m_psingle = 0.001; //TODO

    MagnitudeLogger(const Options &input);

    void log(size_t nexcit, defs::ham_t helem, defs::prob_t prob);

    void synchronize();

};


#endif //M7_MAGNITUDELOGGER_H
