//
// Created by RJA on 20/11/2020.
//

#ifndef M7_FERMIBOSSTOCHASTICPROPAGATOR_H
#define M7_FERMIBOSSTOCHASTICPROPAGATOR_H

#include "StochasticPropagator.h"
#include "src/core/excitgen/BosonCouplingSamplers.h"

#if 0
struct FermiBosStochasticPropagator : public StochasticPropagator {

    BosonCouplingSamplers m_boson_exgen;

    FermiBosStochasticPropagator(const FermiBosHamiltonian& ham, Options& opts):
    StochasticPropagator(ham, opts),
    m_boson_exgen(ham.bc(), opts.nboson_max, m_prng){}

    void diagonal(Wavefunction &m_wf, const size_t &irow) override {
        StochasticPropagator::diagonal(m_wf, irow);
    }

    void off_diagonal(Wavefunction &m_wf, const size_t &irow) override {
        StochasticPropagator::off_diagonal(m_wf, irow);
    }

};


#endif //M7_FERMIBOSSTOCHASTICPROPAGATOR_H
#endif //M7_FERMIBOSSTOCHASTICPROPAGATOR_H
