//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/pchb/HeatBathSamplers.h>
#include "Propagator.h"
#include "FciqmcCalculation.h"


class StochasticPropagator : public Propagator {
    const double &m_min_spawn_mag;
public:
    PRNG m_prng;
    std::unique_ptr<ExcitationGenerator> m_exgen = nullptr;

    StochasticPropagator(FciqmcCalculation *fciqmc) :
            Propagator(fciqmc), m_min_spawn_mag(m_input.min_spawn_mag),
            m_prng(m_input.prng_seed, m_input.prng_ngen) {
        m_exgen = std::unique_ptr<ExcitationGenerator>(new HeatBathSamplers(m_ham.get(), m_prng));
    }

    template<typename T>
    size_t get_nattempt(const T& weight){
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        /*
         * We want to make nattempt = ceil(|weight|) spawning attempts.
         * can't rely on std::abs to provide the right answer in the case of complex arithmetic with
         * a real-valued Hamiltonian and an integral weight, since the sqrt function is not the exact
         * inverse of squaring in finite precision arithmetic!
         */
        if (weight<0) return (-weight) < 1 ? 1 : std::round(-weight);
        else return (weight) < 1 ? 1 : std::round(weight);
    }

    template<typename T>
    size_t get_nattempt(const std::complex<T>& weight) {
        if (m_ham->complex_valued()) return get_nattempt(std::abs(weight));
        else return get_nattempt(consts::real(weight));
    }

    void off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override;

    void diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                  bool flag_deterministic,
                  defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) override;
};

#endif //M7_STOCHASTICPROPAGATOR_H
