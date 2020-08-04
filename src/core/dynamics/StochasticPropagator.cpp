//
// Created by Robert John Anderson on 2020-04-11.
//

#include "StochasticPropagator.h"

void StochasticPropagator::off_diagonal(const DeterminantElement &src_det, const NumericElement<defs::ham_t> &weight,
                                        SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) {
    ASSERT(!consts::float_is_zero(*weight));
    ASSERT(consts::imag(*weight)==0.0 || m_ham->complex_valued())
    m_occ.update(src_det);
    m_vac.update(src_det);
    size_t nattempt = get_nattempt(*weight);
#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "spawn attempts: " << nattempt << std::endl;
#endif
    defs::prob_t prob;
    defs::ham_t helem;
    bool valid = false;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        size_t nexcit = 2 - m_prng.stochastic_round(m_magnitude_logger.m_psingle, 1);
        switch (nexcit) {
            case 1:
                valid = m_exgen->draw_single(src_det, m_dst_det, m_occ, m_vac, prob, helem, m_aconn);
                if (!valid) break;
                ASSERT(prob>=0.0 && prob<=1.0)
                prob*=m_magnitude_logger.m_psingle;
                ASSERT(!consts::float_nearly_zero(prob, 1e-14));
                break;
            case 2:
                // TODO: don't need m_vac for doubles.
                valid = m_exgen->draw_double(src_det, m_dst_det, m_occ, prob, helem, m_aconn);
                if (!valid) break;
                ASSERT(prob>=0.0 && prob<=1.0)
                prob*= 1.0-m_magnitude_logger.m_psingle;
                break;
            default:
                throw std::runtime_error("invalid excitation rank");
        }

#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "EXCITATION GENERATED" << std::endl;
        std::cout << consts::verb << "excitation rank:         " << nexcit << std::endl;
        std::cout << consts::verb << "is valid:                " << string_utils::yn(valid) << std::endl;
#endif

        if (!valid) continue;
        ASSERT(!consts::float_is_zero(prob))
        auto delta = -(*weight / (defs::ham_comp_t) nattempt) * tau() * helem /prob;
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << "probability:             " << prob << std::endl;
        std::cout << consts::verb << "H matrix element:        " << helem << std::endl;
        std::cout << consts::verb << "continuous delta:        " << delta << std::endl;
#endif
        delta = m_prng.stochastic_threshold(delta, m_min_spawn_mag);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << "delta post-thresh:       " << delta << std::endl;
#endif

        ASSERT(consts::floats_equal(delta, -(*weight / (defs::ham_comp_t) nattempt) * tau() * helem /prob)
               || consts::float_is_zero(delta) || consts::float_is_zero(delta-m_min_spawn_mag))

        if (consts::float_is_zero(delta)) continue;
        ASSERT(m_dst_det.nsetbit()==src_det.nsetbit())

        spawn(spawn_list, m_dst_det, delta, flag_initiator, flag_deterministic);
        m_magnitude_logger.log(nexcit, helem, prob);
    }
}

void StochasticPropagator::diagonal(const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight,
                                    bool flag_deterministic, defs::ham_comp_t &delta_square_norm,
                                    defs::ham_comp_t &delta_nw) {

    delta_square_norm -= std::pow(std::abs(*weight), 2);
    delta_nw -= std::abs(*weight);

    if (1||flag_deterministic){
        weight *= 1 - (*hdiag - m_shift) * tau();
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "DETERMINISTIC DEATH" << std::endl;
        std::cout << consts::verb << "new weight:     " << *weight << std::endl;
#endif
    }
    else {
        // the probability that each unit walker will die
        auto death_rate = (*hdiag - m_shift) * tau();
        ASSERT(std::abs(death_rate) < 1)
        weight = m_prng.stochastic_round(*weight, 1.0) * (1 - death_rate);
#ifdef VERBOSE_DEBUGGING
        std::cout << consts::verb << consts::chevs << "STOCHASTIC DEATH" << std::endl;
        std::cout << consts::verb << "death rate:      " << death_rate << std::endl;
        std::cout << consts::verb << "new weight:      " << *weight << std::endl;
        std::cout << consts::verb << "all died:        " << string_utils::yn(consts::float_is_zero(*weight)) << std::endl;
#endif
    }

    delta_square_norm += std::pow(std::abs(*weight), 2);
    delta_nw += std::abs(*weight);
}
