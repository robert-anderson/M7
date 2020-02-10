//
// Created by Robert John Anderson on 2020-01-21.
//

#ifndef M7_COMBINATIONENUMERATOR_H
#define M7_COMBINATIONENUMERATOR_H


#include <cstddef>
#include <string>
#include "../defs.h"
#include "Enumerator.h"

class CombinationEnumerator : public Enumerator<defs::inds> {
    const size_t m_n, m_r;
    std::string m_starting_bitmask, m_bitmask;
    bool m_allfound = false;
public:
    explicit CombinationEnumerator(size_t n, size_t r, Enumerator *subsequent = nullptr);
    inline bool next_element(defs::inds &result) override;
};

#endif //M7_COMBINATIONENUMERATOR_H
