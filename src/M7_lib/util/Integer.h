//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_INTEGER_H
#define M7_UTIL_INTEGER_H

#include "M7_lib/defs.h"

namespace integer {

    /**
     * std::min is not constexpr in C++11 so implement it here
     */
    static constexpr size_t min(size_t i, size_t j) {
        return i < j ? i : j;
    }

    /**
     * std::max is not constexpr in C++11 so implement it here
     */
    static constexpr size_t max(size_t i, size_t j) {
        return i > j ? i : j;
    }

    template<typename T>
    T divceil(T num, T denom) {
        static_assert(std::is_integral<T>::value, "only applicable to integral types");
        return num % denom ? num / denom + 1 : num / denom;
    }

    template<typename T>
    T round_up(T num, T modulo) {
        static_assert(std::is_integral<T>::value, "only applicable to integral types");
        return divceil(num, modulo) * modulo;
    }

    size_t rectmap(size_t irow, size_t icol, size_t ncol);

    void inv_rectmap(size_t &irow, size_t &icol, size_t ncol, size_t flat);

    size_t trigmap(size_t i, size_t j);

    size_t npair(size_t ndim);

    void inv_trigmap(size_t &i, size_t &j, size_t n);

    size_t strigmap(size_t i, size_t j);

    void inv_strigmap(size_t &i, size_t &j, size_t n);

    size_t nspair(size_t ndim);

    size_t factorial(size_t n);

    size_t combinatorial(size_t n, size_t r);

    size_t combinatorial_with_repetition(size_t n, size_t r);

    template<size_t n>
    static typename std::enable_if<n == 0ul, size_t>::type ntup_num(size_t extent) {
        return 1ul;
    }

    template<size_t n>
    static typename std::enable_if<n != 0ul, size_t>::type ntup_num(size_t extent) {
        return extent * ntup_num<n - 1>(extent - 1);
    }

    template<size_t n>
    static size_t ntup(size_t extent) {
        return ntup_num<n>(extent) / ntup_num<n>(n);
    }
}

#endif //M7_UTIL_INTEGER_H
