//
// Created by Robert John Anderson on 2020-02-16.
//

#include <assert.h>
#include <src/utils.h>
#include "VectorCombinationEnumerator.h"

bool VectorCombinationEnumerator::next_element(defs::inds &result) {
    auto tmp = CombinationEnumerator::next_element(result);
    if (!tmp) return false;
    for (auto i{0ul}; i<result.size(); ++i) {
        assert(result[i]<m_vector.size());
        result[i] = m_vector[result[i]];
    }
    return tmp;
}