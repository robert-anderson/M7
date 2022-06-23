//
// Created by Robert John Anderson on 2020-01-17.
//

#ifndef M7_CONSTS_H
#define M7_CONSTS_H

#include <complex>
#include <array>
#include <limits>
#include <typeinfo>

namespace datatype {

    /**
     * generic case: non-complex
     */
    template<typename T>
    struct is_complex_base_t : public std::false_type {};

    /**
     * complex case: is complex
     */
    template<typename T>
    struct is_complex_base_t<std::complex<T>> : public std::true_type {};

    template<typename T>
    struct remove_const_ref_t {
        typedef typename std::remove_reference<typename std::remove_const<T>::type>::type type;
    };

    /**
     * @return
     *  compile time constant bool: true if the template arg is a complex type
     */
    template<typename T>
    constexpr bool is_complex() {
        return is_complex_base_t<typename remove_const_ref_t<T>::type>::value;
    }

    /**
    * the std definition of "arithmetic" type only applies to primitive types, so here the definition is broadened to
    * include the complex numbers.
    */
    template<typename T>
    constexpr bool is_arithmetic() { return is_complex<T>() || std::is_arithmetic<T>::value; }


    template<typename T>
    static std::string name() { return typeid(T).name(); }

    template<>
    std::string name<char>() { return "char"; }

    template<>
    std::string name<short int>() { return "short int"; }

    template<>
    std::string name<int>() { return "int"; }

    template<>
    std::string name<long int>() { return "long int"; }

    template<>
    std::string name<long long int>() { return "long long int"; }

    template<>
    std::string name<unsigned char>() { return "unsigned char"; }

    template<>
    std::string name<unsigned short int>() { return "unsigned short int"; }

    template<>
    std::string name<unsigned int>() { return "unsigned short int"; }

    template<>
    std::string name<unsigned long int>() { return "unsigned long int"; }

    template<>
    std::string name<unsigned long long int>() { return "unsigned long long int"; }

    template<>
    std::string name<float>() { return "float"; }

    template<>
    std::string name<double>() { return "double"; }

    template<>
    std::string name<long double>() { return "long double"; }

    template<>
    std::string name<std::complex<float>>() { return "complex float"; }

    template<>
    std::string name<std::complex<double>>() { return "complex double"; }

    template<>
    std::string name<std::complex<long double>>() { return "complex long double"; }

    template<>
    std::string name<bool>() { return "bool"; }
}

/**
 * defines overloaded and templated functions to deal with complex/real arithmetic appropriately
 */
namespace arith {
    using namespace datatype;

    /**
     * In the generic case, the type of a single complex component "comp" (real or imag part) is simply the type itself
     */
    template<typename T>
    struct comp_base_t {
        typedef T type;
    };

    /**
     * In the complex case, the component type is extracted
     */
    template<typename T>
    struct comp_base_t<std::complex<T>> {
        typedef T type;
    };

    /**
     * convenient definition to extract the component type for all value categories
     */
    template<typename T>
    using comp_t = typename comp_base_t<typename remove_const_ref_t<T>::type>::type;

    /**
     * in the generic case, the complex conjugate is just the input
     */
    template<typename T>
    static constexpr T conj(const T &v) {
        static_assert(!is_complex<T>(), "Complex values should be conjugated by overloads.");
        return v;
    }
    /**
     * in the complex case, the returned value has the imag part negated
     */
    template<typename T>
    static constexpr std::complex<T> conj(const std::complex<T> &v) {
        return std::conj(v);
    }
    /**
     * generic case: real part is just the input
     */
    template<typename T>
    static constexpr T real(const T &v) {
        static_assert(!is_complex<T>(), "Real part of complex values should be taken by overloads.");
        return v;
    }
    /**
     * complex case: real part of number is returned
     */
    template<typename T>
    static constexpr T real(const std::complex<T> &v) {
        return v.real();
    }
    /**
     * generic case: imag part is zero
     */
    template<typename T>
    static constexpr T imag(const T &v) {
        static_assert(!is_complex<T>(), "Imaginary part of complex values should be taken by overloads.");
        return T(0);
    }
    /**
     * complex case: imag part of number is returned
     */
    template<typename T>
    static constexpr T imag(const std::complex<T>& v) {
        return v.imag();
    }
    /**
     * return a const ref to the real part of the complex number
     */
    template<typename T>
    static constexpr const T &real_ref(const std::complex<T> &v) {
        return reinterpret_cast<const std::array<T, 2> &>(v)[0];
    }
    /**
     * return a non-const ref to the real part of the complex number
     */
    template<typename T>
    static constexpr T &real_ref(std::complex<T> &v) {
        return reinterpret_cast<std::array<T, 2> &>(v)[0];
    }
    /**
     * return a const ref to the imag part of the complex number
     */
    template<typename T>
    static constexpr const T &imag_ref(const std::complex<T> &v) {
        return reinterpret_cast<const std::array<T, 2> &>(v)[1];
    }
    /**
     * return a non-const ref to the imag part of the complex number
     */
    template<typename T>
    static constexpr T &imag_ref(std::complex<T> &v) {
        return reinterpret_cast<std::array<T, 2> &>(v)[1];
    }

    /**
     * @tparam T
     *  component type of the complex number
     * @param arg
     *  argument in radians
     * @return
     *  complex number on the unit circle
     */
    template<typename T>
    static std::complex<T> unit_complex(comp_t<T> arg) {
        return {std::cos(arg), std::sin(arg)};
    }

    /**
     * combine real and imag parts into a single complex vector
     * @tparam T
     *  component type
     * @param real
     *  real values
     * @param imag
     *  imag values
     * @param v
     *  complex combination
     */
    template<typename T>
    static void zip(const std::vector<T>& real, const std::vector<T>& imag, std::vector<std::complex<T>>& v){
        auto n = std::min(real.size(), imag.size());
        v.clear();
        v.reserve(n);
        for (size_t i=0ul; i<n; ++i) v.push_back({real[i], imag[i]});
    }

    template<typename T>
    static std::vector<std::complex<T>> zip(const std::vector<T>& real, const std::vector<T>& imag){
        std::vector<std::complex<T>> v;
        zip(real, imag, v);
        return v;
    }
}


/**
 * methods relating to the (approximate) equality of floating point and complex floating point values
 * for floating point comparison we take the same approach as numpy's isclose function
 * https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
 * absolute(a - b) <= (atol + rtol * absolute(b))
 */
namespace fptol {
    using namespace arith;
    /**
     * @tparam T
     *  type of number being compared
     * @return
     *  default relative tolerance for nearly_equal comparisons
     */
    template<typename T>
    static constexpr comp_t<T> default_rtol_near() {return 1e-5;}

    template<typename T>
    static constexpr comp_t<T> default_rtol_near(const T&) {return default_rtol_near<T>();}

    /**
     * @tparam T
     *  type of number being compared
     * @return
     *  default absolute tolerance for nearly_equal comparisons
     */
    template<typename T>
    static constexpr comp_t<T> default_atol_near() {return 1e-8;}

    template<typename T>
    static constexpr comp_t<T> default_atol_near(const T&) {return default_atol_near<T>();}

    template<typename T>
    static constexpr T default_atol_num(const float& v) {return T(0);}
    static constexpr float default_atol_num(const float& v) {return 1e-7;}
    static constexpr double default_atol_num(const double & v) {return 1e-14;}

    template<typename T>
    static constexpr comp_t<T> default_atol_num(const T&) {return default_atol_num<T>();}
    /**
     * a different default atol is used for numeric equality, where one is only leaving room for discrepancies due to
     * non-associativity of floating point ops - useful e.g. in MPI reductions: where the order of operations is
     * non-deterministic, and the number of terms is runtime specifiable
     *
     * this atol is also to be used when comparing against zero, where non-zero rtol is not suitable, although in this
     * case the user can specify another tolerance when the above reasons are not expected to be those responsible
     * discrepancies via nearly_equal
     * @tparam T
     *  type of number being compared
     * @return
     *  default absolute tolerance for numeric_zero comparisons
     */
    template<typename T>
    static constexpr comp_t<T> default_atol_num() {return default_atol_num(comp_t<T>());}

    template<typename T>
    static constexpr bool nearly_equal(T a, T b, T rtol, T atol){
        return std::abs(a-b) <= (atol+rtol*std::abs(b));
    }

    template<typename T>
    static constexpr bool nearly_equal(T a, T b){
        return nearly_equal(a, b, default_rtol_near<T>(), default_atol_near<T>());
    }

    template<typename T>
    static constexpr bool nearly_equal(std::complex<T> a, std::complex<T> b, T rtol, T atol){
        return nearly_equal(a.real(), b.real(), rtol, atol) && nearly_equal(a.real(), b.real(), rtol, atol);
    }

    template<typename T>
    static constexpr bool nearly_equal(std::complex<T> a, std::complex<T> b){
        return nearly_equal(a, b, default_rtol_near<T>(), default_atol_near<T>());
    }

    template<typename T>
    static constexpr bool nearly_zero(T b, comp_t<T> atol){
        return nearly_equal(T(0), b, T(0), atol);
    }

    template<typename T>
    static constexpr bool nearly_real(T b, comp_t<T> atol){
        return numeric_zero(imag(b), atol);
    }

    template<typename T>
    static constexpr bool numeric_equal(T a, T b){
        return nearly_equal(a, b, arith::comp_t<T>(0.0), default_atol_num<T>());
    }

    template<typename T>
    static constexpr bool numeric_zero(T b){
        return nearly_zero(b, default_atol_num<T>());
    }

    template<typename T>
    static constexpr bool numeric_real(T b){
        return numeric_zero(imag(b));
    }

    template<typename T>
    bool numeric_integer(const T &v) {
        static_assert(std::is_floating_point<T>::value, "T must be floating point");
        return numeric_equal(std::round(v), v);
    }
}


#endif //M7_CONSTS_H
