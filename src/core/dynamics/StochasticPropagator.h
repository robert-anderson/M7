//
// Created by Robert John Anderson on 2020-04-11.
//

#ifndef M7_STOCHASTICPROPAGATOR_H
#define M7_STOCHASTICPROPAGATOR_H

#include <src/core/sample/PRNG.h>
#include <src/core/heatbath/HeatBathSampler.h>
#include "Propagator.h"


class StochasticPropagator : public Propagator {
    PRNG m_prng;
    HeatBathSampler m_precomputed_heat_bath_sampler;
public:
    StochasticPropagator(const InputOptions &input, const std::unique_ptr<Hamiltonian> &ham,
                         const RankAllocator<DeterminantElement> &rankAllocator);

public:
    void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                      SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override;

};

#endif //M7_STOCHASTICPROPAGATOR_H
