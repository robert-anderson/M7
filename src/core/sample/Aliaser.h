//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_ALIASER_H
#define M7_ALIASER_H

#include<vector>
#include <src/core/util/defs.h>
#include <stack>
#include <iostream>
#include <src/core/util/utils.h>
#include <src/core/parallel/SharedArray.h>
#include "PRNG.h"

class Aliaser {
    const size_t m_nprob;
    SharedArray<defs::prob_t> m_prob_table;
    SharedArray<size_t> m_alias_table;
    defs::prob_t m_norm;

public:
    Aliaser(const size_t nprob) :
        m_nprob(nprob),
        m_prob_table(m_nprob),
        m_alias_table(m_nprob){
    }

    void update(const defs::prob_t *probs, const size_t nprob) {
        m_norm = std::accumulate(probs, probs+nprob, 0.0);
        std::stack<size_t> smaller;
        std::stack<size_t> larger;
        for (size_t iprob = 0ul; iprob < m_nprob; ++iprob) {
            m_prob_table.set(iprob, probs[iprob] * m_nprob);
            if (m_prob_table[iprob] < m_norm) smaller.push(iprob);
            else larger.push(iprob);
        }
        size_t small, large;
        while (!smaller.empty() && !larger.empty()) {
            small = smaller.top();
            smaller.pop();
            large = larger.top();
            larger.pop();
            m_alias_table.set(small, large);
            m_prob_table.set(large, m_prob_table[large] - (m_norm - m_prob_table[small]));
            if (m_prob_table[large] < m_norm) smaller.push(large);
            else larger.push(large);
        }
    }

    void update(const std::vector<defs::prob_t> &probs) {
        update(probs.data(), probs.size());
    }

    Aliaser(const defs::prob_t *probs, const size_t nprob) :
        Aliaser(nprob) {
        update(probs, nprob);
    }

    Aliaser(const std::vector<defs::prob_t> &probs, const size_t nprob, PRNG &prng) :
        Aliaser(probs.data(), nprob) {}

    Aliaser(const std::vector<defs::prob_t> &probs) :
        Aliaser(probs.data(), probs.size()) {}


    size_t draw(PRNG& prng) const {
        size_t iprob = std::floor(prng.draw_float() * m_nprob);
        ASSERT(iprob < m_nprob);
        if (prng.draw_float() * m_norm < m_prob_table[iprob]) return iprob;
        else return m_alias_table[iprob];
    }

    const defs::prob_t& norm() const {
        return m_norm;
    }

    const size_t& nprob() const {
        return m_nprob;
    }

    defs::prob_t prob(const size_t &i) const {
        return m_prob_table[i]/m_norm;
    }

};


#endif //M7_ALIASER_H
