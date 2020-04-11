//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include "Propagator.h"

class ExactPropagator : public Propagator {

public:
    ExactPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                    const RankAllocator<DeterminantElement> &rank_allocator);

    void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                              SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override ;
};

#endif //M7_EXACTPROPAGATOR_H
