//
// Created by Robert John Anderson on 2020-02-16.
//

#ifndef M7_VECTORCOMBINATIONENUMERATOR_H
#define M7_VECTORCOMBINATIONENUMERATOR_H

#include "CombinationEnumerator.h"
#include "assert.h"

class VectorCombinationEnumerator : public CombinationEnumerator {
    const defs::inds &m_vector;
public:
    VectorCombinationEnumerator(const defs::inds &vector, size_t n, size_t r, Enumerator *subsequent = nullptr):
        CombinationEnumerator(n, r, subsequent), m_vector(vector){}
    VectorCombinationEnumerator(const defs::inds &vector, size_t r, Enumerator *subsequent = nullptr):
        CombinationEnumerator(vector.size(), r, subsequent), m_vector(vector){}
    bool next_element(defs::inds &result) override;
};


#endif //M7_VECTORCOMBINATIONENUMERATOR_H
