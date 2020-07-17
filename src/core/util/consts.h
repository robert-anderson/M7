//
// Created by Robert John Anderson on 2020-01-17.
//

#include <complex>
#include <limits>

#ifndef M7_CONSTS_H
#define M7_CONSTS_H

namespace consts {


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
    constexpr bool is_arithmetic() { return is_complex_t<T>::value || std::is_arithmetic<T>::value; }

    template<typename T>
    constexpr bool is_complex() { return is_complex_t<T>::value; }

    template<typename T>
    static constexpr T conj(const T &v) {
        static_assert(!is_complex<T>(), "Complex values should be conjugated by overloads.");
	   	return v;
    }

    template<typename T>
    static constexpr std::complex<T> conj(const std::complex<T> &v) {
        return std::conj(v);
    }

    template<typename T>
    static constexpr T real(const T &v) {
        static_assert(!is_complex<T>(), "Real part of complex values should be taken by overloads.");
        return v;
    }

    template<typename T>
    static constexpr T real(std::complex<T> v) {
        return std::real(v);
    }

    template<typename T>
    static constexpr T imag(T &v) {
        static_assert(!is_complex<T>(), "Imaginary part of complex values should be taken by overloads.");
        return T(0);
    }

    template<typename T>
    static constexpr T imag(std::complex<T> v) {
        return std::imag(v);
    }

    template<typename T>
    static constexpr T first_quadrant(T &v) {
        return abs(v);
    }

    template<typename T>
    static constexpr std::complex<T> first_quadrant(std::complex<T> v) {
        return std::complex<T>(abs(v.real()), abs(v.imag()));
    }

    template<typename T>
    static constexpr T l1_norm(T &v) {
        return abs(v);
    }

    template<typename T>
    static constexpr T l1_norm(std::complex<T> v) {
        return std::abs(v.real())+std::abs(v.imag());
    }

    template<typename T>
    static constexpr T component_norm(T &v) {
        return std::abs(v);
    }

    template<typename T>
    static constexpr std::complex<T> component_norm(std::complex<T> v) {
        return std::complex<T>(std::abs(v.real()), std::abs(v.imag()));
    }

    template<typename T>
    static T real_log(T &v) {
        return std::log(std::abs(v));
    }

    template<typename T>
    static std::complex<T> real_log(std::complex<T> v) {
        return std::complex<T>(std::log(std::abs(v.real())), std::log(std::abs(v.imag())));
    }

    template<typename T>
    static constexpr T real_ratio(T &v1, T &v2) {
        return v1/v2;
    }

    template<typename T>
    static constexpr std::complex<T> real_ratio(std::complex<T> v1, std::complex<T> v2) {
        return std::complex<T>(v1.real()/v2.real(), v1.imag()/v2.imag());
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
        return std::abs(v)<eps;
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
    static constexpr bool floats_nearly_equal(T v1, T v2, typename component_t<T>::type eps=1e-12){
        return float_nearly_zero(v1-v2, eps);
    }

    const double pi = std::atan(1.0)*4;
    const double two_pi = 2*pi;
    const double sqrt2 = std::sqrt(2.0);
    const double invsqrt2 = 1.0/sqrt2;

}


#endif //M7_CONSTS_H
