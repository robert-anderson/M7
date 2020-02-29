//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_PRNG_H
#define M7_PRNG_H

#include <random>
#include <assert.h>

class PRNG {
    std::mt19937 m_mt19937;
public:
    PRNG(const size_t &seed):m_mt19937(seed){
        assert(m_mt19937.min()==0u);
        assert(m_mt19937.max()==~0u);
    }

    uint32_t draw_uint(){
        return m_mt19937();
    }

    double draw_float(){
        return double(draw_uint())/(1+std::mt19937::max());
    }

};


#endif //M7_PRNG_H
