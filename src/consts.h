//
// Created by Robert John Anderson on 2020-01-17.
//

#include <complex>

#ifndef M7_CONSTS_H
#define M7_CONSTS_H

namespace consts {

    template<typename T>
    static constexpr bool is_floating_point = std::is_floating_point<T>::value;

    template<typename T>
    struct is_complex_t : public std::false_type {};

    template<typename T>
    struct is_complex_t<std::complex<T>> : public std::true_type {};

    template<typename T>
    struct is_complex_t<const std::complex<T>> : public std::true_type {};

    template<typename T>
    struct is_complex_t<std::complex<T>&> : public std::true_type {};

    template<typename T>
    struct is_complex_t<const std::complex<T>&> : public std::true_type {};

    template<typename T>
    constexpr bool is_complex() { return is_complex_t<T>::value; }

    template<typename T>
    static constexpr T conj(T &v) {
        if constexpr (is_complex<T>()) return std::conj(v);
        else return v;
    }

    template<typename T>
    static constexpr T real(T &v) {
        return v;
    }

    template<typename T>
    static constexpr T real(std::complex<T> v) {
        return std::real(v);
    }

    template<typename T>
    static constexpr T imag(T &v) {
        return T(0);
    }

    template<typename T>
    static constexpr T imag(std::complex<T> v) {
        return std::imag(v);
    }

    template<typename T>
    struct component_t {
        typedef T type;
    };

    template<typename T>
    struct component_t<std::complex<T>> {
        typedef T type;
    };

    template<typename T>
    struct component_t<const std::complex<T>&> {
        typedef T type;
    };

    template<typename T>
    static constexpr bool float_nearly_zero(T v, typename component_t<T>::type eps){
        return abs(v)<eps;
    }

    template<typename T>
    static constexpr bool float_is_zero(T v){
        return float_nearly_zero(v, std::numeric_limits<typename component_t<T>::type>::epsilon());
    }

    template<typename T>
    static constexpr bool floats_equal(T v1, T v2){
        return float_is_zero(v1-v2);
    }

    template<typename T>
    static constexpr bool floats_nearly_equal(T v1, T v2, typename component_t<T>::type eps=1e-14){
        return float_nearly_zero(v1-v2, eps);
    }

}

#endif //M7_CONSTS_H