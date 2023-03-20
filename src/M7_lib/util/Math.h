//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_MATH_H
#define M7_UTIL_MATH_H

#include <complex>
#include <numeric>
#include "M7_lib/defs.h"


namespace math {

    const double pi = std::atan(1.0) * 4;
    const double two_pi = 2 * pi;
    const double sqrt2 = std::sqrt(2.0);
    const double invsqrt2 = 1.0 / sqrt2;

    /**
     * exact raising to integer power by recursive squaring
     * @tparam T
     *  type of base and result
     * @param x
     *  number being exponentiated
     * @param exp
     *  integer exponent
     * @return
     *  x to the power y
     */
    template<typename T>
    T pow(T x, uint_t exp) {
        if (!exp) return 1;
        if (exp == 1) return x;

        auto tmp = pow(x, exp / 2);
        if (exp & 1ul) return x * tmp * tmp;
        return tmp * tmp;
    }

    /**
     * exponentiation in the case that the exponent is a compile time constant integer
     * @tparam exp
     *  integer exponent
     * @tparam T
     *  type of base and result
     * @param x
     *  number being exponentiated
     * @return
     *  x to the power exp
     *
     */
    template<uint_t exp, typename T=void>
    static typename std::enable_if<exp == 0ul, T>::type pow(T /*x*/) {
        return 1ul;
    }

    template<uint_t exp, typename T=void>
    static typename std::enable_if<exp != 0ul, T>::type pow(T x) {
        return x * pow<exp - 1, T>(x);
    }

    template<typename T>
    T l1_norm(const T* ptr, uint_t n) {
        return std::accumulate(ptr, ptr+n, T{}, [](const T& tot, const T& v){
            return tot+std::abs(v);
        });
    }

    template<typename T>
    T l2_norm(const T* ptr, uint_t n) {
        return std::sqrt(std::accumulate(ptr, ptr+n, T{}, [](const T& tot, const T& v){
            return tot+pow<2>(std::abs(v));
        }));
    }

    /**
     * behaves like std::clamp from C++17 but returns copy instead of reference
     */
    template<typename T>
    T clamp(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v ? hi : v);
    }

    template<typename T>
    T phase(const T& v) {
        return v/std::abs(v);
    }
}

#endif //M7_UTIL_MATH_H
