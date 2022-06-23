//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_ALIASER_H
#define M7_ALIASER_H

#include <vector>
#include <M7_lib/defs.h>
#include <stack>
#include <iostream>
#include <M7_lib/parallel/SharedMatrix.h>
#include <M7_lib/parallel/MPIAssert.h>
#include "PRNG.h"

struct Aliaser {
    const size_t m_nrow, m_nprob;
private:
    SharedMatrix<defs::prob_t> m_prob_table;
    SharedMatrix<size_t> m_alias_table;
    SharedArray<defs::prob_t> m_norm;

public:
    Aliaser(size_t nrow, size_t nprob);

    void update(size_t irow, const defs::prob_t *probs, size_t nprob);

    void update(size_t irow, const std::vector<defs::prob_t> &probs);

    size_t draw(size_t irow, PRNG &prng) const;

    defs::prob_t norm(size_t irow) const;
};

class SingleAliaser : public Aliaser {
public:
    SingleAliaser(size_t nprob);

    SingleAliaser(const std::vector<defs::prob_t> &probs);

    void update(const std::vector<defs::prob_t> &probs);

    size_t draw(PRNG &prng) const;

    defs::prob_t norm() const;
};

#endif //M7_ALIASER_H
