//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <M7_lib/sample/PRNG.h>
#include <M7_lib/excitgen/ExcitGen.h>
#include <M7_lib/excitgen/UniformSingles.h>
#include <M7_lib/excitgen/HeatBathDoubles.h>
#include <M7_lib/excitgen/ExcitGenGroup.h>
#include <M7_lib/excitgen2/ExcitGenGroup2.h>

#include "Propagator.h"

class StochasticPropagator : public Propagator {
protected:
    PRNG m_prng;
    ExcitGenGroup m_excit_gens;
    ExcitGenGroup2 m_excit_gen_group;
    MagnitudeLogger m_mag_log;
    const double &m_min_spawn_mag;

    template<typename T>
    static T phase(const T& weight) {
        return weight < 0.0 ? -1.0: 1.0;
    }

    template<typename T>
    static std::complex<T> phase(const std::complex<T>& weight) {
        return weight/std::abs(weight);
    }

public:
    StochasticPropagator(const Hamiltonian &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    void diagonal(Wavefunction &m_wf, const size_t& ipart) override;

    template<typename T>
    size_t get_nattempt(const T &weight) {
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        return m_prng.stochastic_round(std::abs(weight), 1.0);
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
