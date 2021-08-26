//
// Created by Robert John Anderson on 2020-02-20.
//

#include "Aliaser.h"

Aliaser::Aliaser(const size_t &nrow, const size_t &nprob) :
        m_nrow(nrow),
        m_nprob(nprob),
        m_prob_table(nrow, m_nprob),
        m_alias_table(nrow, m_nprob),
        m_norm(nrow) {
    DEBUG_ONLY(m_nrow);
}

Aliaser::Aliaser(const std::vector<defs::prob_t> &probs) : Aliaser(1, probs.size()) {
    update(0, probs);
}

void Aliaser::update(const size_t &irow, const defs::prob_t *probs, const size_t nprob) {
    ASSERT(irow < m_nrow)
    m_norm.set(irow, std::accumulate(probs, probs + nprob, 0.0));
    std::stack<size_t> smaller;
    std::stack<size_t> larger;
    for (size_t iprob = 0ul; iprob < m_nprob; ++iprob) {
        m_prob_table.set(irow, iprob, probs[iprob] * m_nprob);
        if (m_prob_table.get(irow, iprob) < m_norm[irow]) smaller.push(iprob);
        else larger.push(iprob);
    }
    size_t small, large;
    while (!smaller.empty() && !larger.empty()) {
        small = smaller.top();
        smaller.pop();
        large = larger.top();
        larger.pop();
        m_alias_table.set(irow, small, large);
        m_prob_table.set(irow, large,
                         m_prob_table.get(irow, large) - (m_norm[irow] - m_prob_table.get(irow, small)));
        if (m_prob_table.get(irow, large) < m_norm[irow]) smaller.push(large);
        else larger.push(large);
    }
}

void Aliaser::update(const size_t &irow, const std::vector<defs::prob_t> &probs) {
    update(irow, probs.data(), probs.size());
}

void Aliaser::update(const std::vector<defs::prob_t> &probs) {
    ASSERT(m_nrow==1)
    update(0, probs);
}

size_t Aliaser::draw(const size_t &irow, PRNG &prng) const {
    ASSERT(irow < m_nrow)
    size_t iprob = prng.draw_uint(m_nprob);
    ASSERT(iprob < m_nprob)
    if (prng.draw_float() * m_norm[irow] < m_prob_table.get(irow, iprob)) return iprob;
    else return m_alias_table.get(irow, iprob);
}

size_t Aliaser::draw(PRNG &prng) const {
    ASSERT(m_nrow==1)
    return draw(0, prng);
}

const defs::prob_t &Aliaser::norm(const size_t &irow) const {
    ASSERT(irow < m_nrow)
    return m_norm[irow];
}

const size_t &Aliaser::nprob() const {
    return m_nprob;
}

defs::prob_t Aliaser::prob(const size_t &irow, const size_t &iprob) const {
    return m_prob_table.get(irow, iprob) / m_norm[irow];
}
