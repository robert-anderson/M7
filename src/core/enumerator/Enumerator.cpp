//
// Created by rja on 21/03/2021.
//

#include "Enumerator.h"

enums::Enumerator::Enumerator(size_t n, size_t r, size_t nv) : m_n(n), m_r(r), m_nv(nv), m_v(r, 0ul) {}

const size_t &enums::Enumerator::operator[](const size_t &i) const {
    return m_v[i];
}

enums::CombinationsWithRepetition::CombinationsWithRepetition(size_t n, size_t r) : Enumerator(n, r, integer_utils::combinatorial(n, r)+n) {
    m_v.back() = ~0ul;
}

bool enums::CombinationsWithRepetition::next() {
    m_v.back() += 1;
    for (size_t i = m_v.size() - 1; i != ~0ul; --i) {
        if (m_v[i] >= m_n) {
            if (!i) {
                m_v.assign(m_v.size(), 0);
                m_v.back() = ~0ul;
                return false;
            }
            m_v[i - 1] += 1;
            for (size_t k = 0ul; k <= i; ++k) {
                m_v[k] = m_v[i - 1];
            }
        }
    }
    return true;
}

enums::CombinationsDistinct::CombinationsDistinct(size_t n, size_t r) : Enumerator(n, r, integer_utils::combinatorial(n, r)), m_starting_bitmask(r, 1){
    ASSERT(m_n >= m_r);
    m_starting_bitmask.resize(n, 0);
    m_bitmask = m_starting_bitmask;
}

bool enums::CombinationsDistinct::next() {
    if (m_allfound) {
        m_bitmask = m_starting_bitmask;
        m_allfound = false;
        return false;
    }
    size_t ipos = 0ul;
    for (size_t i=0ul; i<m_n; ++i) {
        if (m_bitmask[i]){
            m_v[ipos] = i;
            if (ipos++==m_r) break;
        }
    }
    m_allfound = !std::prev_permutation(m_bitmask.begin(), m_bitmask.end());
#ifndef NDEBUG
    auto is_ordered = [&](){
        for (size_t i=1ul; i<m_v.size(); ++i) if (m_v[i]<=m_v[i-1]){return false;}
        return true;
    };
    ASSERT(is_ordered());
#endif
    return true;
}

enums::PermutationsWithRepetition::PermutationsWithRepetition(size_t n, size_t r) : Enumerator(n, r, std::pow(n, r)) {
    m_v.back() = ~0ul;
}

bool enums::PermutationsWithRepetition::next() {
    // find index to increment;
    for (size_t i = m_r - 1; i != ~0ul; --i) {
        if (m_v[i] + 1 < m_n) {
            m_v[i]++;
            // reset higher positions
            std::fill(m_v.begin() + i +1, m_v.end(), 0ul);
            return true;
        }
    }
    std::fill(m_v.begin(), m_v.end(), 0);
    m_v.back() = ~0ul;
    return false;
}
