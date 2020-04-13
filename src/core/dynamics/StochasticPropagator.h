//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/heatbath/HeatBathSampler.h>
#include "Propagator.h"
#include "FciqmcCalculation.h"


class StochasticPropagator : public Propagator {
public:
    PRNG m_prng;
    HeatBathSampler m_precomputed_heat_bath_sampler;
    StochasticPropagator(FciqmcCalculation *fciqmc) :
        Propagator(fciqmc),
        m_prng(fciqmc->m_input.prng_seed, fciqmc->m_input.prng_ngen),
        m_precomputed_heat_bath_sampler(fciqmc->m_ham.get(), m_prng) {}

    void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override {

        defs::ham_comp_t largest_spawned_magnitude = 0.0;
        assert(!consts::float_is_zero(*weight));

        size_t nattempt = std::ceil(std::abs(*weight));
        auto &det_sampler = m_precomputed_heat_bath_sampler.det_sampler->get(0);
        det_sampler.update(determinant);

        for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
            det_sampler.draw();
            if (det_sampler.single_generated()) {
                assert(!consts::float_is_zero(det_sampler.get_single_prob()));
                auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau *
                             (m_ham->get_element_1(det_sampler.get_single()) / det_sampler.get_single_prob());
                if (consts::float_is_zero(delta)) continue;
                spawn(spawn_list, det_sampler.get_single_dst_det(), delta, largest_spawned_magnitude, flag_initiator);
            }
            if (det_sampler.double_generated()) {
                assert(!consts::float_is_zero(det_sampler.get_double_prob()));
                auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau *
                             (m_ham->get_element_2(det_sampler.get_double()) / det_sampler.get_double_prob());
                if (consts::float_is_zero(delta)) continue;
                spawn(spawn_list, det_sampler.get_double_dst_det(), delta, largest_spawned_magnitude, flag_initiator);
            }
        }
    }
};

#endif //M7_STOCHASTICPROPAGATOR_H
