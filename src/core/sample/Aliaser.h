//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_ALIASER_H
#define M7_ALIASER_H

#include <vector>
#include <src/defs.h>
#include <stack>
#include <iostream>
#include <src/utils.h>
#include <src/core/thread/PrivateStore.h>
#include "PRNG.h"

class Aliaser {
    const size_t m_nprob;
    PrivateStore<PRNG> &m_prng;
    std::vector<defs::prob_t> m_prob_table;
    defs::inds m_alias_table;
    defs::prob_t m_norm;

public:
    Aliaser(const size_t nprob, PrivateStore<PRNG> &prng) :
        m_nprob(nprob),
        m_prng(prng),
        m_prob_table(m_nprob, 0.0),
        m_alias_table(m_nprob, 0ul){}

    void update(const defs::prob_t *probs, const size_t nprob) {
        m_norm = std::accumulate(probs, probs+nprob, 0.0);
        std::stack<size_t> smaller;
        std::stack<size_t> larger;
        for (size_t iprob = 0ul; iprob < m_nprob; ++iprob) {
            m_prob_table[iprob] = probs[iprob] * m_nprob;
            if (m_prob_table[iprob] < m_norm) smaller.push(iprob);
            else larger.push(iprob);
        }
        size_t small, large;
        while (!smaller.empty() && !larger.empty()) {
            small = smaller.top();
            smaller.pop();
            large = larger.top();
            larger.pop();
            m_alias_table[small] = large;
            m_prob_table[large] -= (m_norm - m_prob_table[small]);
            if (m_prob_table[large] < m_norm) smaller.push(large);
            else larger.push(large);
        }
    }

    void update(const std::vector<defs::prob_t> &probs) {
        update(probs.data(), probs.size());
    }

    Aliaser(const defs::prob_t *probs, const size_t nprob, PrivateStore<PRNG> &prng) :
        Aliaser(nprob, prng) {
        update(probs, nprob);
    }

    Aliaser(const std::vector<defs::prob_t> &probs, const size_t nprob, PrivateStore<PRNG> &prng) :
        Aliaser(probs.data(), nprob, prng) {}

    Aliaser(const std::vector<defs::prob_t> &probs, PrivateStore<PRNG> &prng) :
        Aliaser(probs.data(), probs.size(), prng) {}


    size_t draw() const {
        size_t iprob = std::floor(m_prng.get(0).draw_float() * m_nprob);
        ASSERT(iprob >= 0);
        ASSERT(iprob < m_nprob);
        if (m_prng.get(0).draw_float() * m_norm < m_prob_table[iprob]) return iprob;
        else return m_alias_table[iprob];
    }

    defs::prob_t norm() const {
        return m_norm;
    }

};


#endif //M7_ALIASER_H
