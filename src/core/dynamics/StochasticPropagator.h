//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/excitgen/ExcitationGenerator.h>
#include <src/core/excitgen/UniformSingles.h>
#include <src/core/excitgen/HeatBathDoubles.h>
#include <src/core/sample/WeightedDrawer.h>
#include <src/core/excitgen/BosonExcitationGenerator.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {

    void add_boson_excitgen(const Hamiltonian &ham);

protected:
    PRNG m_prng;
    std::vector<std::unique_ptr<ExcitationGenerator>> m_exgens;
    const double &m_min_spawn_mag;
    std::unique_ptr<WeightedDrawer> m_exgen_drawer;

public:
    StochasticPropagator(const Hamiltonian &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt);

    void diagonal(Wavefunction &m_wf, const size_t& ipart) override;

    template<typename T>
    size_t get_nattempt(const T &weight) {
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");

#ifndef ENABLE_CEILING_SPAWN_ATTEMPTS
        return m_prng.stochastic_round(std::abs(weight), 1.0);
#else
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
        if (m_ham.complex_valued()) return get_nattempt(std::abs(weight));
        else return get_nattempt(consts::real(weight));
    }

    void off_diagonal(Wavefunction &wf, const size_t& ipart) override;

    bool is_exact() const override;
};

#endif //M7_STOCHASTICPROPAGATOR_H
