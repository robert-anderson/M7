//
// Created by rja on 27/02/2020.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H


#include <src/sample/PRNG.h>
#include "Propagator.h"

class StochasticPropagator : public Propagator {
    PRNG m_prng;
public:
    StochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                         const RankAllocator<Determinant> &rank_allocator);


    void off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                      const NumericView<bool> flag_deterministic, const NumericView<bool> flag_initiator,
                      TableArray<SpawnList> &spawn_list) const override;
};


#endif //M7_STOCHASTICPROPAGATOR_H
