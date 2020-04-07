//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_PRNG_H
#define M7_PRNG_H

#include <random>
#include <assert.h>
#include <algorithm>

class PRNG {
    std::mt19937 m_mt19937;
    std::vector<uint32_t> m_data;
    std::vector<uint32_t>::const_iterator m_it = m_data.end();
public:
    PRNG(const size_t &seed, const size_t &block_size):
    m_mt19937(seed), m_data(block_size) {
        assert(m_mt19937.min()==0u);
        assert(m_mt19937.max()==~0u);
    }

    void refresh() {
        std::uniform_int_distribution<uint32_t> dist(m_mt19937.min(), m_mt19937.max());
        std::generate(m_data.begin(), m_data.end(), [&](){ return dist(m_mt19937); });
        m_it = m_data.begin();
    }

    uint32_t draw_uint(){
        if (m_it==m_data.end()) refresh();
        return *m_it++;
    }

    double draw_float(){
        return double(draw_uint())/(1ul+std::mt19937::max());
    }

};


#endif //M7_PRNG_H
