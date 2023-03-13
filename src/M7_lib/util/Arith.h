//
// Created by rja on 23/06/22.
//

#ifndef M7_ARITH_H
#define M7_ARITH_H

#include <array>
#include <complex>
#include <vector>
#include "M7_lib/defs.h"
#include "Datatype.h"


/**
 * defines overloaded and templated functions to deal with complex/real arithmetic appropriately
 */
namespace arith {

    /**
     * real or complex number with a given component type
     */
    template <typename comp_t, bool real>
    using num_t = typename std::conditional<real, comp_t, std::complex<comp_t>>::type;

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
    using comp_t = typename comp_base_t<typename dtype::remove_const_ref_t<T>::type>::type;

    /**
     * in the generic case, the complex conjugate is just the input
     */
    template<typename T>
    static constexpr T conj(const T &v) {
        static_assert(dtype::is_arithmetic<T>(), "only conjugation of arithmetic or complex types makes sense");
        static_assert(!dtype::is_complex<T>(), "Complex values should be conjugated by overloads.");
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
        static_assert(dtype::is_arithmetic<T>(), "only real part of arithmetic or complex types makes sense");
        static_assert(!dtype::is_complex<T>(), "Real part of complex values should be taken by overloads.");
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
    static constexpr T imag(const T&) {
        static_assert(dtype::is_arithmetic<T>(), "only imag part of arithmetic or complex types makes sense");
        static_assert(!dtype::is_complex<T>(), "Imaginary part of complex values should be taken by overloads.");
        return T(0);
    }

    /**
     * complex case: imag part of number is returned
     */
    template<typename T>
    static constexpr T imag(const std::complex<T> &v) {
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
    static void zip(const v_t<T> &real, const v_t<T> &imag, v_t<std::complex<T>> &v) {
        auto n = std::min(real.size(), imag.size());
        v.clear();
        v.reserve(n);
        for (uint_t i = 0ul; i < n; ++i) v.push_back({real[i], imag[i]});
    }

    template<typename T>
    static v_t<std::complex<T>> zip(const v_t<T> &real, const v_t<T> &imag) {
        v_t<std::complex<T>> v;
        zip(real, imag, v);
        return v;
    }
}



#endif //M7_ARITH_H
