//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_MATH_H
#define M7_UTIL_MATH_H

#include <complex>
#include <numeric>
#include <random>
#include "M7_lib/defs.h"
#include "Arith.h"


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

    /**
     * log of geometric mean of the absolute value of given data
     */
    template<typename T>
    T log_geo_mean(const T* v, uint_t n) {
        T tot = 0.0;
        for (auto ptr=v; ptr!=v+n; ++ptr) tot+=std::log(std::abs(*ptr));
        return tot / n;
    }

    template<typename T>
    T log_geo_mean(const v_t<T>& v) {
        return log_geo_mean(v.data(), v.size());
    }

    /**
     * geometric mean of the absolute value of given data
     * mean(v) = (prod_i |v_i|) / sqrt(n)
     *         = exp(sum_i log v_i / n)
     */
    template<typename T>
    T geo_mean(const T* v, uint_t n) {
        return std::exp(log_geo_mean(v, n));
    }

    template<typename T>
    T geo_mean(const v_t<T>& v) {
        return geo_mean(v.data(), v.size());
    }

    template <typename T>
    class Logarithm {
        T m_magnitude;
        bool m_sign;
    public:
        Logarithm(const T& magnitude, bool sign): m_magnitude(magnitude), m_sign(!dtype::is_complex<T>() && sign) {}
        explicit Logarithm(const T& v): Logarithm(std::log(std::abs(v)), arith::real(v) < 0) {}

        Logarithm& operator += (const Logarithm<T>& other) {
            m_magnitude += other.m_magnitude;
            m_sign ^= other.m_sign;
            return *this;
        }

        Logarithm operator + (const Logarithm<T>& other) const {
            auto tmp = *this;
            tmp += other;
            return tmp;
        }

        Logarithm& operator -= (const Logarithm<T>& other) {
            m_magnitude -= other.m_magnitude;
            m_sign ^= other.m_sign;
            return *this;
        }

        Logarithm operator - (const Logarithm<T>& other) const {
            auto tmp = *this;
            tmp -= other;
            return tmp;
        }

        Logarithm& operator *= (const T& v) {
            *this += Logarithm<T>(v);
            return *this;
        }

        Logarithm operator * (const T& v) const {
            auto tmp = *this;
            tmp *= v;
            return tmp;
        }

        Logarithm& operator /= (const T& v) {
            *this -= Logarithm<T>(v);
            return *this;
        }

        Logarithm operator / (const T& v) const {
            auto tmp = *this;
            tmp /= v;
            return tmp;
        }

        T exp() const {
            return m_sign ? -std::exp(m_magnitude) : std::exp(m_magnitude);
        }
    };


    template<typename T>
    Logarithm<T> logarithm(const T& magnitude, bool sign) {return {magnitude, sign};}
    template<typename T>
    Logarithm<T> logarithm(const T& v) {return Logarithm<T>(v);}
}

#endif //M7_UTIL_MATH_H
