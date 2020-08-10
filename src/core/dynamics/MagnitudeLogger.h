//
// Created by rja on 18/05/2020.
//

#ifndef M7_MAGNITUDELOGGER_H
#define M7_MAGNITUDELOGGER_H


#include<cstddef>
#include <src/core/util/defs.h>
#include <src/core/parallel/Reducible.h>
#include <src/core/io/Options.h>
#include <src/core/parallel/Epoch.h>

class MagnitudeLogger {
    const Options &m_input;
    size_t m_nsingle = 0;
    size_t m_ndouble = 0;
    // highest magnitudes
    Reducible<defs::ham_comp_t> m_hi_mag_single;
    Reducible<defs::ham_comp_t> m_hi_mag_double;

    Epoch m_enough_singles_for_dynamic_tau;
    Epoch m_enough_doubles_for_dynamic_tau;

public:
    // the recommended timestep based on the hi_mag and the maximum acceptable bloom
    defs::prob_t m_psingle;
    double m_tau;

    MagnitudeLogger(const Options &input, defs::prob_t m_psingle);

    MagnitudeLogger(const Options &input, size_t nsite, size_t nelec):
    MagnitudeLogger(input,1/(integer_utils::combinatorial(2*nsite-nelec, 2)/double(2*nsite-nelec)+1))
    {
        std::cout << "Magnitude Logger initialized" << std::endl;
    }

    void log(size_t nexcit, defs::ham_t helem, defs::prob_t prob);

    void synchronize(size_t icycle);

};


#endif //M7_MAGNITUDELOGGER_H
