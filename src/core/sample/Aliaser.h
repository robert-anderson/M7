//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_ALIASER_H
#define M7_ALIASER_H

#include <vector>
#include <src/defs.h>
#include <stack>
#include <iostream>
#include <src/core/util/utils.h>
#include <src/core/parallel/SharedMatrix.h>
#include <src/core/parallel/MPIAssert.h>
#include "PRNG.h"

class Aliaser {
    const size_t m_nrow, m_nprob;
    SharedMatrix<defs::prob_t> m_prob_table;
    SharedMatrix<size_t> m_alias_table;
    SharedArray<defs::prob_t> m_norm;

public:
    Aliaser(const size_t &nrow, const size_t &nprob);

    Aliaser(const size_t &nprob) : Aliaser(1, nprob) {}

    Aliaser(const std::vector<defs::prob_t> &probs);

    void update(const size_t &irow, const defs::prob_t *probs, const size_t nprob);

    void update(const size_t &irow, const std::vector<defs::prob_t> &probs);

    void update(const std::vector<defs::prob_t> &probs);

    size_t draw(const size_t &irow, PRNG &prng) const;

    size_t draw(PRNG &prng) const;

    const defs::prob_t &norm(const size_t &irow) const;

    const size_t &nprob() const;
};

#endif //M7_ALIASER_H
