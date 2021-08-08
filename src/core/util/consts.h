//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_CONSTS_H
#define M7_CONSTS_H

#include <complex>
#include <array>
#include <limits>
#include <typeinfo>

namespace consts {

    template<typename T>
    struct is_complex_t : public std::false_type {
    };

    template<typename T>
    struct is_complex_t<std::complex<T>> : public std::true_type {
    };

    template<typename T>
    struct is_complex_t<const std::complex<T>> : public std::true_type {
    };

    template<typename T>
    struct is_complex_t<std::complex<T> &> : public std::true_type {
    };

    template<typename T>
    struct is_complex_t<const std::complex<T> &> : public std::true_type {
    };


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
    static constexpr T real(const std::complex<T> &v) {
        return v.real();
    }

    template<typename T>
    static constexpr T imag(const T &v) {
        static_assert(!is_complex<T>(), "Imaginary part of complex values should be taken by overloads.");
        return T(0);
    }

    template<typename T>
    static constexpr T imag(const std::complex<T>& v) {
        return v.imag();
    }

    template<typename T>
    static constexpr const T &real_ref(const std::complex<T> &v) {
        return reinterpret_cast<const std::array<T, 2> &>(v)[0];
    }

    template<typename T>
    static constexpr T &real_ref(std::complex<T> &v) {
        return reinterpret_cast<std::array<T, 2> &>(v)[0];
    }

    template<typename T>
    static constexpr const T &imag_ref(const std::complex<T> &v) {
        return reinterpret_cast<const std::array<T, 2> &>(v)[1];
    }

    template<typename T>
    static constexpr T &imag_ref(std::complex<T> &v) {
        return reinterpret_cast<std::array<T, 2> &>(v)[1];
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
        return std::abs(v.real()) + std::abs(v.imag());
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
        return v1 / v2;
    }

    template<typename T>
    static constexpr std::complex<T> real_ratio(std::complex<T> v1, std::complex<T> v2) {
        return std::complex<T>(v1.real() / v2.real(), v1.imag() / v2.imag());
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
    struct component_t<const std::complex<T> &> {
        typedef T type;
    };

    template<typename T>
    using real_t = typename component_t<T>::type;

    template<typename T>
    static constexpr bool float_nearly_zero(T v, typename component_t<T>::type eps) {
        return std::abs(v) < eps;
    }

    template<typename T>
    static constexpr bool float_is_zero(T v) {
        return float_nearly_zero(v, std::numeric_limits<typename component_t<T>::type>::epsilon());
    }

    template<typename T>
    static constexpr bool floats_equal(T v1, T v2) {
        return float_is_zero(v1 - v2);
    }

    template<typename T>
    static constexpr bool floats_nearly_equal(T v1, T v2, typename component_t<T>::type eps = 1e-12) {
        return float_nearly_zero(v1 - v2, eps);
    }

    const double pi = std::atan(1.0) * 4;
    const double two_pi = 2 * pi;
    const double sqrt2 = std::sqrt(2.0);
    const double invsqrt2 = 1.0 / sqrt2;

    // for debugging output
    const std::string verb = "\t[VERBOSE]  ";
    const std::string chevs = " >>> ";

    template<typename T>
    static std::string type_name() { return typeid(T).name(); }

    template<>
    std::string type_name<char>() { return "char"; }

    template<>
    std::string type_name<short int>() { return "short int"; }

    template<>
    std::string type_name<int>() { return "int"; }

    template<>
    std::string type_name<long int>() { return "long int"; }

    template<>
    std::string type_name<long long int>() { return "long long int"; }

    template<>
    std::string type_name<unsigned char>() { return "unsigned char"; }

    template<>
    std::string type_name<unsigned short int>() { return "unsigned short int"; }

    template<>
    std::string type_name<unsigned int>() { return "unsigned short int"; }

    template<>
    std::string type_name<unsigned long int>() { return "unsigned long int"; }

    template<>
    std::string type_name<unsigned long long int>() { return "unsigned long long int"; }

    template<>
    std::string type_name<float>() { return "float"; }

    template<>
    std::string type_name<double>() { return "double"; }

    template<>
    std::string type_name<long double>() { return "long double"; }

    template<>
    std::string type_name<std::complex<float>>() { return "complex float"; }

    template<>
    std::string type_name<std::complex<double>>() { return "complex double"; }

    template<>
    std::string type_name<std::complex<long double>>() { return "complex long double"; }

    template<>
    std::string type_name<bool>() { return "bool"; }

}


#endif //M7_CONSTS_H
