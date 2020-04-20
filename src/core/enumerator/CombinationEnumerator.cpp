//
// Created by Robert John Anderson on 2020-01-21.
//

#include "CombinationEnumerator.h"
#include <algorithm>
#include <iostream>
#include "src/defs.h"

CombinationEnumerator::CombinationEnumerator(size_t n, size_t r, Enumerator* subsequent) :
Enumerator(subsequent), m_n(n), m_r(r), m_starting_bitmask(r, 1)
{
    ASSERT(m_n >= m_r);
    m_starting_bitmask.resize(n, 0);
    m_bitmask = m_starting_bitmask;
}

bool CombinationEnumerator::next_element(defs::inds &result) {
    ASSERT(result.size() == m_r);
    if (m_allfound) {
        m_bitmask = m_starting_bitmask;
        m_allfound = false;
        return false;
    }
    size_t ipos = 0ul;
    for (size_t i=0ul; i<m_n; ++i) {
        if (m_bitmask[i]){
            result[ipos] = i;
            if (ipos++==m_r) break;
        }
    }
    auto is_ordered = [](const defs::inds inds){
        for (size_t i=1ul; i<inds.size(); ++i) if (inds[i]<=inds[i-1]){return false;}
        return true;
    };
    m_allfound = !std::prev_permutation(m_bitmask.begin(), m_bitmask.end());
    ASSERT(is_ordered(result));
    return true;
}
