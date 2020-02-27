//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include "DeterministicPropagator.h"

class ExactPropagator : public DeterministicPropagator {

public:
    ExactPropagator(const std::unique_ptr<Hamiltonian> &ham, const RankAllocator<Determinant> &rankAllocator,
                    double tau, defs::ham_comp_t shift);

    void off_diagonal(const Determinant &determinant, const NumericView<defs::ham_t> &weight,
                      const NumericView<bool> flag_deterministic, const NumericView<bool> flag_initiator,
                      TableArray<SpawnList> &spawn_list) const override;
};


#endif //M7_EXACTPROPAGATOR_H
