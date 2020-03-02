//
// Created by rja on 27/02/2020.
//

#include <src/heatbath/DeterminantSampler.h>
#include "StochasticPropagator.h"
#include "src/consts.h"


StochasticPropagator::StochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                           const RankAllocator<Determinant> &rank_allocator) :
        Propagator(input, ham, rank_allocator), m_prng(input.prng_seed),
        m_precomputed_heat_bath_sampler(*ham) {}

void StochasticPropagator::off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                                        const NumericView<bool> flag_deterministic,
                                        const NumericView<bool> flag_initiator,
                                        TableArray<SpawnList> &spawn_list) {
    size_t nattempt = std::ceil(std::abs(*weight));
    DeterminantSampler determinant_sampler(m_precomputed_heat_bath_sampler, determinant);
    for (size_t iattempt=0ul; iattempt<nattempt; ++iattempt){
        auto excit = determinant_sampler.draw(m_prng);
        if (!excit.m_single.is_null()) {
            assert(!consts::float_is_zero(excit.m_single.m_helement));
            assert(!consts::float_is_zero(excit.m_single.m_prob));
            auto excited = excit.m_single.get_connection();
            auto delta = -(*weight / (defs::ham_comp_t) nattempt) * m_tau *
                         (m_ham->get_element(excited, determinant) / excit.m_single.m_prob);
            add_to_spawn_list(excited, delta, flag_initiator, spawn_list);
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
