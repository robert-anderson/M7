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

    void off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override {
        ASSERT(!consts::float_is_zero(*weight));
        m_occ.update(src_det);
        m_vac.update(src_det);
        size_t nattempt = std::ceil(std::abs(*weight));
        defs::prob_t prob;
        defs::ham_t helem;
        bool valid = false;

        for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
            size_t nexcit = 2 - m_prng.stochastic_round(m_magnitude_logger.m_psingle, 1);
            switch (nexcit) {
                case 1:
                    valid = m_exgen->draw_single(src_det, m_dst_det, m_occ, m_vac, prob, helem, m_aconn);
                    prob*=m_magnitude_logger.m_psingle;
                    ASSERT(!consts::float_nearly_zero(prob, 1e-14));
                    break;
                case 2:
                    // TODO: don't need m_vac for doubles.
                    valid = m_exgen->draw_double(src_det, m_dst_det, m_occ, prob, helem, m_aconn);
                    prob*= 1.0-m_magnitude_logger.m_psingle;
                    break;
            }
            if (!valid) continue;

            ASSERT(!consts::float_is_zero(prob))
            auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau * helem /prob;
            delta = m_prng.stochastic_threshold(delta, m_min_spawn_mag);
            if (consts::float_is_zero(delta)) continue;
            ASSERT(m_dst_det.nsetbit()==src_det.nsetbit())
            spawn(spawn_list, m_dst_det, delta, flag_initiator, flag_deterministic);
            m_magnitude_logger.log(nexcit, helem, prob);
        }
    }

    void diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                  bool flag_deterministic,
                  defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) override {

        delta_square_norm -= std::pow(std::abs(*weight), 2);
        delta_nw -= std::abs(*weight);

        if (flag_deterministic){
            weight *= 1 - (*hdiag - m_shift) * m_tau;
        }
        else {
            // the probability that each unit walker will die
            auto death_rate = (*hdiag - m_shift) * m_tau;
            ASSERT(std::abs(death_rate) < 1)
            weight = m_prng.stochastic_round(*weight, 1.0) * (1 - death_rate);
        }

        delta_square_norm += std::pow(std::abs(*weight), 2);
        delta_nw += std::abs(*weight);
    }
};

#endif //M7_STOCHASTICPROPAGATOR_H
