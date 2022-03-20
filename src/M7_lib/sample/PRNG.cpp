//
// Created by Robert John Anderson on 2020-02-20.
//

#include <M7_lib/parallel/MPIWrapper.h>
#include <M7_lib/parallel/MPIAssert.h>
#include "PRNG.h"

PRNG::PRNG(const size_t &seed, const size_t &block_size) :
        m_data(block_size, 0u), m_seed(seed+mpi::irank()) {
    ASSERT(block_size > 0);
    m_i = m_data.size();
}

void PRNG::refresh() {
    std::mt19937 mt19937(m_seed+(mpi::irank()+1)*m_nrefresh);
    std::uniform_int_distribution<uint32_t> dist(mt19937.min(), mt19937.max());
    std::generate(m_data.begin(), m_data.end(), [&]() { return dist(mt19937); });
    m_i = 0;
    m_nrefresh++;
}

uint32_t PRNG::draw_uint() {
    if (m_i == m_data.size()) refresh();
    return m_data[m_i++];
}

uint32_t PRNG::draw_uint(uint32_t modular_base) {
    return draw_uint() % modular_base;
}

uint32_t PRNG::draw_uint(uint32_t min, uint32_t max){
    DEBUG_ASSERT_LT(min, max, "invalid limits");
    return min + draw_uint(max-min);
}

double PRNG::draw_float() {
    return double(draw_uint()) / (1ul + std::mt19937::max());
}

uint32_t PRNG::draw_uint_linear_bias(uint32_t modular_base) {
    return std::floor(modular_base*std::sqrt(draw_float()));
}
