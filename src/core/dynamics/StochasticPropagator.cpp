//
// Created by Robert John Anderson on 2020-04-11.
//

#include "StochasticPropagator.h"

#if 0
void
StochasticPropagator::off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                                   SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) {

    size_t nattempt = std::ceil(std::abs(*weight));
    DeterminantSampler determinant_sampler(m_precomputed_heat_bath_sampler, determinant, m_prng);
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        determinant_sampler.draw();
        if (determinant_sampler.single_generated()) {
            assert(!consts::float_is_zero(determinant_sampler..m_single.m_helement));
            assert(!consts::float_is_zero(excit.m_single.m_prob));
            auto excited = excit.m_single.get_connection();
            auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau *
                         (m_ham->get_element(excited, determinant) / excit.m_single.m_prob);
            add_to_spawn_list(excited, delta, flag_initiator, spawn_list);
        }
    }
    if (!excit.m_double.is_null()) {
        assert(!consts::float_is_zero(excit.m_double.m_helement));
        assert(!consts::float_is_zero(excit.m_double.m_prob));
        auto excited = excit.m_double.get_connection();
        auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau *
                     (m_ham->get_element(excited, determinant) / excit.m_double.m_prob);
        add_to_spawn_list(excited, delta, flag_initiator, spawn_list);
    }
}

}

StochasticPropagator::StochasticPropagator(
    const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
    const RankAllocator<DeterminantElement> &rankAllocator) :
    Propagator(input, ham, rankAllocator), m_prng(input.prng_seed, 1000),
    m_precomputed_heat_bath_sampler(*ham.get()) {}

#endif