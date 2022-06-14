//
// Created by Robert J. Anderson on 17/08/2020.
//

#ifndef M7_KRAMERSSECTOROCCUPATION_H
#define M7_KRAMERSSECTOROCCUPATION_H

#if 0

#include <M7_lib/parallel/Reducible.h>
#include <M7_lib/util/defs.h>
#include <M7_lib/basis/DeterminantField.h>

class KramersSectorOccupation {
    const size_t m_nelec;
    std::vector<Reducible<defs::wf_comp_t>> m_sum;

public:
    explicit KramersSectorOccupation(size_t nelec);

    explicit KramersSectorOccupation(const DeterminantElement& det);

    ~KramersSectorOccupation();

    void add(const DeterminantElement& det, const defs::wf_t& weight);

};


#endif //M7_KRAMERSSECTOROCCUPATION_H
#endif //M7_KRAMERSSECTOROCCUPATION_H
