//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/excitgen/ExcitGen.h>
#include <src/core/excitgen/UniformSingles.h>
#include <src/core/excitgen/HeatBathDoubles.h>
#include <src/core/excitgen/ExcitGenGroup.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {
protected:
    PRNG m_prng;
    ExcitGenGroup m_excit_gens;
    MagnitudeLogger m_mag_log;
    const double &m_min_spawn_mag;

public:
    StochasticPropagator(const Hamiltonian &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    void diagonal(Wavefunction &m_wf, const size_t& ipart) override;

    template<typename T>
    size_t get_nattempt(const T &weight) {
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        return m_prng.stochastic_round(std::abs(weight), 1.0);
#if 0
        /*
         * We want to make nattempt = ceil(|weight|) spawning attempts.
         * can't rely on std::abs to provide the right answer in the case of complex arithmetic with
         * a real-valued FermionHamiltonian and an integral weight, since the sqrt function is not the exact
         * inverse of squaring in finite precision arithmetic!
         */
        if (weight < 0) return (-weight) < 1 ? 1 : std::round(-weight);
        else return (weight) < 1 ? 1 : std::round(weight);
#endif
    }

    template<typename T>
    size_t get_nattempt(const std::complex<T> &weight) {
        return get_nattempt(std::abs(weight));
    }

    void off_diagonal(Wavefunction &wf, const size_t& ipart) override;

    size_t nexcit_gen() const override;

    std::vector<defs::prob_t> exlvl_probs() const override;

    void update(const size_t &icycle, const Wavefunction &wf) override;

};

#endif //M7_STOCHASTICPROPAGATOR_H
