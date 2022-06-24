//
// Created by Robert John Anderson on 2020-02-20.
//

#ifndef M7_PRNG_H
#define M7_PRNG_H

#include <random>
#include <algorithm>
#include <complex>
#include <M7_lib/defs.h>
#include <M7_lib/parallel/MPIAssert.h>

class PRNG {
    typedef std::vector<uint32_t> U;
    U m_data;
    uint_t m_i;
    const uint_t m_seed;
    uint_t m_nrefresh = 0;
public:
    PRNG(const uint_t &seed, const uint_t &block_size);

    void refresh();

    uint32_t draw_uint();

    uint32_t draw_uint(uint32_t modular_base);

    uint32_t draw_uint(uint32_t min, uint32_t max);

    /**
     * sample integers according to a linear distribution (favoring larger integers in the range)
     * @param modular_base
     * @return
     */
    uint32_t draw_uint_linear_bias(uint32_t modular_base);

    double draw_float();

    /**
     * @tparam T
     *  non-complex floating point type to round
     * @param v
     *  value to round
     * @param magnitude
     *  scalar about which to round v
     * @param prob
     *  probability of the lower magnitude outcome
     * @return
     *  the lower magnitude result with probability given by prob, or the higher magnitude result with probability given
     *  by 1-prob
     */
    template<typename T>
    T stochastic_round(const T &v, const T &magnitude, defs::prob_t& prob) {
        static_assert(std::is_floating_point<T>::value, "stochastic round is only applicable to floating point types");
        DEBUG_ASSERT_GE(magnitude, 1e-14, "magnitude should be non-negative and non-zero");
        prob = v / magnitude;
        const long int_ratio = prob;
        if (v >= 0) {
            prob = prob - double(int_ratio);
            return double(draw_float() < prob ? (int_ratio + 1) : int_ratio) * magnitude;
        }
        else {
            prob = double(int_ratio) - prob;
            return (draw_float() < prob ? (int_ratio - 1) : int_ratio) * magnitude;
        }
    }

    template<typename T>
    std::complex<T> stochastic_round(const std::complex<T> &v, const T &magnitude, defs::prob_t& prob) {
        return stochastic_round(std::abs(v), magnitude, prob) * v / std::abs(v);
    }

    /**
     * overload for situations where the probability is not required by the calling scope
     */
    template<typename T>
    T stochastic_round(const T &v, const arith::comp_t<T> &magnitude) {
        defs::prob_t prob;
        return stochastic_round(v, magnitude, prob);
    }

    /**
     * @tparam T
     *  non-complex floating point type to round
     * @param v
     *  value to round
     * @param magnitude
     *  scalar about which to round v
     * @param prob
     *  probability that the result is rounded up or unmodified
     * @return
     *  result of stochastic_round only if |v| is less than magnitude, else return v unmodified with certainty
     */
    template<typename T>
    T stochastic_threshold(const T &v, const arith::comp_t<T> &magnitude, defs::prob_t& prob) {
        if (std::abs(v) < magnitude) return stochastic_round(v, magnitude, prob);
        prob = 1.0;
        return v;
    }

    /**
     * overload for situations where the probability is not required by the calling scope
     */
    template<typename T>
    T stochastic_threshold(const T &v, const arith::comp_t<T> &magnitude) {
        defs::prob_t prob;
        return stochastic_threshold(v, magnitude, prob);
    }

};


#endif //M7_PRNG_H
