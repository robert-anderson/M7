//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#if 0
class StochasticPropagator : public Propagator {

public:
    ExactPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                    const RankAllocator<DeterminantElement> &rank_allocator);

    void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override ;
};


#endif //M7_STOCHASTICPROPAGATOR_H
#endif //M7_STOCHASTICPROPAGATOR_H
