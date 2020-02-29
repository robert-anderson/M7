//
// Created by rja on 27/02/2020.
//

#include "StochasticPropagator.h"


StochasticPropagator::StochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                                           const RankAllocator<Determinant> &rank_allocator) :
        Propagator(input, ham, rank_allocator), m_prng(input.prng_seed) {}

void StochasticPropagator::off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                                        const NumericView<bool> flag_deterministic,
                                        const NumericView<bool> flag_initiator,
                                        TableArray<SpawnList> &spawn_list) const {
    size_t nattempt = std::ceil(std::abs(*weight));
    for (size_t iattempt=0ul; iattempt<nattempt; ++iattempt){

    }
}
