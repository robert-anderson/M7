//
// Created by Robert John Anderson on 2020-01-21.
//

#include "CombinationEnumerator.h"
#include <algorithm>
#include <assert.h>

CombinationEnumerator::CombinationEnumerator(size_t n, size_t r, Enumerator* subsequent) :
Enumerator(subsequent), m_n(n), m_r(r), m_starting_bitmask(r, 1)
{
    m_starting_bitmask.resize(n, 0);
    m_bitmask = m_starting_bitmask;
}

bool CombinationEnumerator::next_element(defs::inds &result) {
    assert(result.size() == m_r);
    if (m_allfound) {
        m_bitmask = m_starting_bitmask;
        m_allfound = false;
        return false;
    }
    auto ipos{0ul};
    for (auto i{0ul}; i<m_n; ++i) {
        if (m_bitmask[i]){
            result[ipos] = i;
            if (ipos++==m_r) break;
        }
    }
    m_allfound = !std::prev_permutation(m_bitmask.begin(), m_bitmask.end());
    return true;
}
