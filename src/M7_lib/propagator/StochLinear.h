//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_LINEAR_STOCH_PROPAGATOR_H
#define M7_LINEAR_STOCH_PROPAGATOR_H

#include <M7_lib/sample/PRNG.h>
#include <M7_lib/excitgen/ExcitGenGroup.h>
#include "M7_lib/propagator/Propagator.h"

class StochLinear : public Propagator {
protected:
    PRNG m_prng;
    ExcitGenGroup m_excit_gen_group;
    MagnitudeLogger m_mag_log;
    const double m_min_spawn_mag;
    const double m_min_death_mag;

    template<typename T>
    static T phase(const T& weight) {
        return weight < 0.0 ? -1.0: 1.0;
    }

    template<typename T>
    static std::complex<T> phase(const std::complex<T>& weight) {
        return weight/std::abs(weight);
    }

public:
    StochLinear(const Hamiltonian &ham, const conf::Document &opts, const wf::Vectors& wf);

    void diagonal(wf::Vectors &wf, Walker& walker, uint_t ipart) override;

    template<typename T>
    uint_t get_nattempt(const T &weight) {
        static_assert(std::is_floating_point<T>::value, "template arg must be floating point");
        return m_prng.stochastic_round(std::abs(weight), 1.0);
    }

    template<typename T>
    uint_t get_nattempt(const std::complex<T> &weight) {
        return get_nattempt(std::abs(weight));
    }

    void off_diagonal(wf::Vectors &wf, const Walker& walker, uint_t ipart, bool initiator) override;

    uint_t ncase_excit_gen() const override;

    v_t<prob_t> excit_gen_case_probs() const override;

    void update(uint_t icycle, const wf::Vectors &wf, const wf::Refs& refs) override;

    const ExcitGenGroup& excit_gen_group() const {
        return m_excit_gen_group;
    }

    hash::digest_t checksum_() const override;

};

#endif //M7_LINEAR_STOCH_PROPAGATOR_H
