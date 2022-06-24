//
// Created by Robert John Anderson on 2020-02-20.
//

#include "Aliaser.h"

Aliaser::Aliaser(uint_t nrow, uint_t nprob) :
        m_nrow(nrow),
        m_nprob(nprob),
        m_prob_table(m_nrow, m_nprob),
        m_alias_table(m_nrow, m_nprob),
        m_norm(m_nrow) {}

void Aliaser::update(uint_t irow, const prob_t *probs, const uint_t nprob) {
    DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
    m_norm.set(irow, std::accumulate(probs, probs + nprob, 0.0));
    std::stack<uint_t> smaller;
    std::stack<uint_t> larger;
    for (uint_t iprob = 0ul; iprob < m_nprob; ++iprob) {
        m_prob_table.set(irow, iprob, probs[iprob] * m_nprob);
        if (m_prob_table.get(irow, iprob) < m_norm[irow]) smaller.push(iprob);
        else larger.push(iprob);
    }
    uint_t small, large;
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

void Aliaser::update(uint_t irow, const std::vector<prob_t> &probs) {
    update(irow, probs.data(), probs.size());
}

uint_t Aliaser::draw(uint_t irow, PRNG &prng) const {
    DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
    uint_t iprob = prng.draw_uint(m_nprob);
    DEBUG_ASSERT_LT(iprob, m_nprob, "prob index OOB");
    if (prng.draw_float() * m_norm[irow] < m_prob_table.get(irow, iprob)) return iprob;
    else return m_alias_table.get(irow, iprob);
}

prob_t Aliaser::norm(uint_t irow) const {
    DEBUG_ASSERT_LT(irow, m_nrow, "row index OOB");
    return m_norm[irow];
}

SingleAliaser::SingleAliaser(uint_t nprob) : Aliaser(1, nprob){}

SingleAliaser::SingleAliaser(const std::vector<prob_t> &probs) : Aliaser(1, probs.size()) {
    update(probs);
}

void SingleAliaser::update(const std::vector<prob_t> &probs) {
    Aliaser::update(0ul, probs);
}

uint_t SingleAliaser::draw(PRNG &prng) const {
    return Aliaser::draw(0ul, prng);
}

prob_t SingleAliaser::norm() const {
    return Aliaser::norm(0ul);
}
