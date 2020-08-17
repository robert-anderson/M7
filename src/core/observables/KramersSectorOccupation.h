//
// Created by rja on 17/08/2020.
//

#ifndef M7_KRAMERSSECTOROCCUPATION_H
#define M7_KRAMERSSECTOROCCUPATION_H

#include <src/core/parallel/Reducible.h>
#include "src/core/util/defs.h"
#include "src/core/table/DeterminantField.h"

class KramersSectorOccupation {
    const size_t m_nelec, m_nsite;
    std::vector<Reducible<defs::wf_comp_t>> m_sum;

public:
    KramersSectorOccupation(size_t nelec, size_t nsite);

    explicit KramersSectorOccupation(const DeterminantElement& det);

    ~KramersSectorOccupation();

    size_t add(const DeterminantElement& det, const defs::wf_t& weight);

};


#endif //M7_KRAMERSSECTOROCCUPATION_H
