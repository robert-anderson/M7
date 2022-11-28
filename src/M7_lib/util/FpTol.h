//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_FPTOL_H
#define M7_FPTOL_H

#include "Arith.h"

/**
 * methods relating to the (approximate) equality of floating point and complex floating point values.
 * for floating point comparison we take the same approach as numpy's isclose function
 * https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
 * absolute(a - b) <= (atol + rtol * absolute(b))
 */
namespace fptol {
    /**
     * @return
     *  default relative tolerance for near_equal comparisons
     */
    template<typename T>
    static constexpr T default_rtol(const T&) {return {};}
    static constexpr float default_rtol(const float&) {return 1e-5;}
    static constexpr double default_rtol(const double&) {return 1e-5;}

    template<typename T>
    static constexpr arith::comp_t<T> default_rtol(const std::complex<T>& z) {
        return default_rtol(arith::real_ref(z));
    }
    /**
     * @return
     *  default absolute tolerance for near_equal comparisons
     */
    template<typename T>
    static constexpr T default_atol(const T&) {return {};}
    static constexpr float default_atol(const float&) {return 1e-8;}
    static constexpr double default_atol(const double&) {return 1e-8;}

    template<typename T>
    static constexpr arith::comp_t<T> default_atol(const std::complex<T>& z) {
        return default_atol(arith::real_ref(z));
    }

    /**
     * @return
     *  default absolute tolerance for near_zero comparisons
     */
    template<typename T>
    static constexpr T default_ztol(const T&) {return {};}
    static constexpr float default_ztol(const float&) {return 1e-12;}
    static constexpr double default_ztol(const double&) {return 1e-12;}

    template<typename T>
    static constexpr arith::comp_t<T> default_ztol(const std::complex<T>& z) {
        return default_ztol(arith::real_ref(z));
    }


    template<typename T>
    static constexpr bool near_equal(T a, T b, T rtol, T atol){
        return std::abs(a-b) <= (atol+rtol*std::abs(b));
    }

    template<typename T>
    static constexpr bool near_equal(T a, T b){
        return near_equal(a, b, default_rtol(a), default_atol(b));
    }

    template<typename T>
    static constexpr bool near_equal(std::complex<T> a, std::complex<T> b, T rtol, T atol){
        return near_equal(a.real(), b.real(), rtol, atol) && near_equal(a.real(), b.real(), rtol, atol);
    }

    template<typename T>
    static constexpr bool near_equal(std::complex<T> a, std::complex<T> b){
        return near_equal(a, b, default_rtol(a), default_atol(a));
    }

    template<typename T>
    static constexpr bool near_zero(T b, arith::comp_t<T> ztol) {
        return near_equal(T{}, b, arith::comp_t<T>{}, ztol);
    }

    template<typename T>
    static constexpr bool near_zero(T b) {
        return near_zero(b, default_ztol(b));
    }

    template<typename T>
    static constexpr bool near_real(T b, arith::comp_t<T> atol) {
        return near_zero(arith::imag(b), atol);
    }

    template<typename T>
    static constexpr bool near_real(T b){
        return near_zero(arith::imag(b));
    }

    template<typename T>
    bool near_integer(const T &v) {
        static_assert(std::is_floating_point<T>::value, "T must be floating point");
        return near_equal(std::round(v), v);
    }

    template<typename T>
    bool near_integer(const std::complex<T> &v) {
        return near_real(v) && near_integer(arith::real_ref(v));
    }
}


#endif //M7_FPTOL_H
