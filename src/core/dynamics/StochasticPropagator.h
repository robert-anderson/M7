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
    PrivateStore<PRNG> m_prng;
    std::unique_ptr<ExcitationGenerator> m_exgen = nullptr;
    PrivateStore<Determinant> m_dst_det;
    PrivateStore<OccupiedOrbitals> m_occ;
    PrivateStore<VacantOrbitals> m_vac;
    PrivateStore<AntisymConnection> m_anticonn;

    StochasticPropagator(FciqmcCalculation *fciqmc) :
            Propagator(fciqmc), m_min_spawn_mag(m_input.min_spawn_mag),
            m_prng(m_input.prng_seed, m_input.prng_ngen),
            m_dst_det(m_ham->nsite()), m_occ(m_dst_det.get()), m_vac(m_dst_det.get()), m_anticonn(m_dst_det.get()) {
        m_exgen = std::unique_ptr<ExcitationGenerator>(new HeatBathSamplers(m_ham.get(), m_prng));
    }

    void off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override {
        ASSERT(!consts::float_is_zero(*weight));
        m_occ.get().update(src_det);
        m_vac.get().update(src_det);
        size_t nattempt = std::ceil(std::abs(*weight));
        defs::prob_t prob;
        defs::ham_t helem;
        bool valid;

        for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
            size_t nexcit = 2 - m_prng.get().stochastic_round(m_magnitude_logger.m_psingle, 1);
            switch (nexcit) {
                case 1:
                    valid = m_exgen->draw_single(src_det, m_dst_det.get(), m_occ.get(), m_vac.get(), prob, helem, m_anticonn.get());
                    prob*=m_magnitude_logger.m_psingle;
                    ASSERT(!consts::float_nearly_zero(prob, 1e-14));
                    break;
                case 2:
                    // TODO: don't need m_vac for doubles.
                    valid = m_exgen->draw_double(src_det, m_dst_det.get(), m_occ.get(), prob, helem, m_anticonn.get());
                    prob*= 1.0-m_magnitude_logger.m_psingle;
                    break;
            }
            if (!valid) continue;

            ASSERT(!consts::float_is_zero(prob))
            auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau * helem /prob;
            delta = m_prng.get().stochastic_threshold(delta, m_min_spawn_mag);
            if (consts::float_is_zero(delta)) continue;
            spawn(spawn_list, m_dst_det.get(), delta, flag_initiator);
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
            weight = m_prng.get().stochastic_round(*weight, 1.0) * (1 - death_rate);
        }

        delta_square_norm += std::pow(std::abs(*weight), 2);
        delta_nw += std::abs(*weight);
    }
};

#endif //M7_STOCHASTICPROPAGATOR_H
