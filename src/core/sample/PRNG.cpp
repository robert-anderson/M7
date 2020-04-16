//
// Created by Robert John Anderson on 2020-02-20.
//

#include "PRNG.h"

PRNG::PRNG(const size_t &seed, const size_t &block_size) :
    m_seed(seed), m_data(block_size, 0u) {
    m_it = m_data.end();
}

void PRNG::refresh() {
    std::mt19937 mt19937(m_seed+m_data.back());
    std::uniform_int_distribution<uint32_t> dist(mt19937.min(), mt19937.max());
    std::generate(m_data.begin(), m_data.end(), [&](){ return dist(mt19937); });
    m_it = m_data.begin();
}

uint32_t PRNG::draw_uint() {
    if (m_it==m_data.end()) refresh();
    return *m_it++;
}

double PRNG::draw_float() {
    return double(draw_uint())/(1ul+std::mt19937::max());
}
