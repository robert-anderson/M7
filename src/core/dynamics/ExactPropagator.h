//
// Created by rja on 27/02/2020.
//

#ifndef M7_EXACTPROPAGATOR_H
#define M7_EXACTPROPAGATOR_H

#include "Propagator.h"

class ExactPropagator : public Propagator {

public:
    ExactPropagator(FciqmcCalculation *fciqmc);

    void off_diagonal(const DeterminantElement &determinant, const NumericElement<defs::ham_t> &weight,
                              SpawnList &spawn_list, bool flag_deterministic, bool flag_initiator) override ;
};

#endif //M7_EXACTPROPAGATOR_H
