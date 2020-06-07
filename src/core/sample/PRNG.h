//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_PRNG_H
#define M7_PRNG_H

#include <random>
#include <algorithm>
#include <complex>
#include <src/core/thread/AlignedAllocator.h>
#include "src/core/util/defs.h"

class PRNG {
    typedef std::vector<uint32_t, AlignedAllocator<uint32_t, defs::cache_line_size>> U;
    U m_data;
    U::iterator m_it;
    const size_t m_seed;
public:
    PRNG(const size_t &seed, const size_t &block_size);

    void refresh();

    uint32_t draw_uint();

    uint32_t draw_uint(uint32_t);

    double draw_float();

    template <typename T>
    T stochastic_round(const T& v, const double& magnitude){
        ASSERT(magnitude > 0);
        static_assert(std::is_floating_point<T>::value, "Stochastic round is only applicable to floating point types");
        const double ratio = v/magnitude;
        const long int_ratio = ratio;
        if (ratio>=0){
            if (draw_float()<(ratio-int_ratio)) return (int_ratio+1)*magnitude;
            else return int_ratio*magnitude;
        }
        else {
            if (draw_float()<(int_ratio-ratio)) return (int_ratio-1)*magnitude;
            else return int_ratio*magnitude;
        }
    }

    template <typename T>
    std::complex<T> stochastic_round(const std::complex<T>& v, const double& magnitude) {
        return stochastic_round(std::abs(v), magnitude)*v/std::abs(v);
    }

    template <typename T>
    T stochastic_threshold(const T& v, const double& magnitude){
        if (std::abs(v)<magnitude) return stochastic_round(v, magnitude);
        return v;
    }

};


#endif //M7_PRNG_H
