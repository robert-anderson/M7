//
// Created by RJA on 20/11/2020.
//

#include "WeightedDrawer.h"

void WeightedDrawer::update_cumprobs() {
    m_cumprobs.front() = m_probs.front();
    for(size_t iprob = 1ul; iprob<m_probs.size(); ++iprob)
        m_cumprobs[iprob] = m_cumprobs[iprob-1]+m_probs[iprob];
    ASSERT(consts::floats_nearly_equal(1.0, m_cumprobs.back()));
}

void WeightedDrawer::set_one(size_t iprob, defs::prob_t prob) {
    ASSERT(iprob+1==m_probs.size() || iprob+2==m_probs.size());
    m_probs[iprob] = prob;
    if (iprob+2==m_probs.size()) {
        /*
         * final value is set to conserve probability
         */
        m_probs.back() = 1.0;
        for (size_t i = 0ul; i<=iprob; ++i) m_probs.back()-=m_probs[i];
        ASSERT(m_probs.back()>=0.0);
    }
    ASSERT(consts::floats_nearly_equal(1.0, std::accumulate(m_probs.cbegin(), m_probs.cend(), 0.0)));
}

WeightedDrawer::WeightedDrawer(size_t nprob, PRNG &prng) : m_probs(nprob, 1.0/nprob), m_cumprobs(nprob), m_prng(prng){
    update_cumprobs();
}

size_t WeightedDrawer::draw() {
    auto tmp = m_prng.draw_float();
    for (size_t iprob=0ul; iprob<m_probs.size(); ++iprob) if (tmp<m_cumprobs[iprob]) return iprob;
    throw std::runtime_error("Probabilities unnormalized, or prng is out of [0.0, 1.0) range");
}

const defs::prob_t &WeightedDrawer::prob(size_t iprob) {
    return m_probs[iprob];
}
