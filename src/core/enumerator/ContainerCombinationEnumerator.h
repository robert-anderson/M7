//
// Created by Robert John Anderson on 2020-02-16.
//

#ifndef M7_CONTAINERCOMBINATIONENUMERATOR_H
#define M7_CONTAINERCOMBINATIONENUMERATOR_H

#include "CombinationEnumerator.h"
#include "assert.h"

template<typename T>
class ContainerCombinationEnumerator : public CombinationEnumerator {
    const T &m_container;
public:
    ContainerCombinationEnumerator(const T &container, size_t n, size_t r, Enumerator *subsequent = nullptr) :
        CombinationEnumerator(n, r, subsequent), m_container(container) {}

    ContainerCombinationEnumerator(const T &container, size_t r, Enumerator *subsequent = nullptr) :
        CombinationEnumerator(container.size(), r, subsequent), m_container(container) {}

    bool next_element(defs::inds &result) override {
        auto tmp = CombinationEnumerator::next_element(result);
        if (!tmp) return false;
        for (size_t i = 0ul; i < result.size(); ++i) {
            assert(result[i] < m_container.size());
            result[i] = m_container[result[i]];
        }
        return tmp;
    }
};


#endif //M7_CONTAINERCOMBINATIONENUMERATOR_H
