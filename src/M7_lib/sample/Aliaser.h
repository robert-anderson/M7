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
    const uint_t m_nrow, m_nprob;
private:
    SharedMatrix<prob_t> m_prob_table;
    SharedMatrix<uint_t> m_alias_table;
    SharedArray<prob_t> m_norm;

public:
    Aliaser(uint_t nrow, uint_t nprob);

    void update_(uint_t irow, const prob_t *probs, uint_t nprob);

    void update_(uint_t irow, const v_t<prob_t> &probs);

    uint_t draw(uint_t irow, PRNG &prng) const;

    prob_t norm(uint_t irow) const;
};

class SingleAliaser : public Aliaser {
public:
    SingleAliaser(uint_t nprob);

    SingleAliaser(const v_t<prob_t> &probs);

    void update_(const v_t<prob_t> &probs);

    uint_t draw(PRNG &prng) const;

    prob_t norm() const;
};

#endif //M7_ALIASER_H
