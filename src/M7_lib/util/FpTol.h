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
     * @tparam T
     *  type of number being compared
     * @return
     *  default relative tolerance for nearly_equal comparisons
     */
    template<typename T>
    static constexpr arith::comp_t<T> default_rtol_near() {return 1e-5;}

    template<typename T>
    static constexpr arith::comp_t<T> default_rtol_near(const T&) {return default_rtol_near<T>();}

    /**
     * @tparam T
     *  type of number being compared
     * @return
     *  default absolute tolerance for nearly_equal comparisons
     */
    template<typename T>
    static constexpr arith::comp_t<T> default_atol_near() {return 1e-8;}

    template<typename T>
    static constexpr arith::comp_t<T> default_atol_near(const T&) {return default_atol_near<T>();}

    template<typename T>
    static constexpr T default_atol_num(const float&) {return T(0);}
    static constexpr float default_atol_num(const float&) {return 1e-7;}
    static constexpr double default_atol_num(const double&) {return 1e-14;}

    template<typename T>
    static constexpr arith::comp_t<T> default_atol_num(const T&) {return default_atol_num<T>();}
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
    static constexpr arith::comp_t<T> default_atol_num() {return default_atol_num(arith::comp_t<T>());}

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
    static constexpr bool nearly_zero(T b, arith::comp_t<T> atol){
        return nearly_equal(T(0), b, T(0), atol);
    }

    template<typename T>
    static constexpr bool nearly_real(T b, arith::comp_t<T> atol){
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
        return numeric_zero(arith::imag(b));
    }

    template<typename T>
    bool numeric_integer(const T &v) {
        static_assert(std::is_floating_point<T>::value, "T must be floating point");
        return numeric_equal(std::round(v), v);
    }
}


#endif //M7_FPTOL_H
