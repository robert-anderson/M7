//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_INTEGER_H
#define M7_UTIL_INTEGER_H

#include "M7_lib/defs.h"

using namespace defs;

namespace integer {

    /**
     * std::min is not constexpr in C++11 so implement it here
     */
    static constexpr uint_t min(uint_t i, uint_t j) {
        return i < j ? i : j;
    }

    /**
     * std::max is not constexpr in C++11 so implement it here
     */
    static constexpr uint_t max(uint_t i, uint_t j) {
        return i > j ? i : j;
    }

    /**
     * the "ceiling" of an integral division. if the remainder is zero, return the quotient, else return quotient + 1
     */
    template<typename T>
    T divceil(T num, T denom) {
        static_assert(std::is_integral<T>::value, "only applicable to integral types");
        const auto quot = num / denom;
        return quot + ((quot*denom)!=num);
    }

    template<typename T>
    T round_up(T num, T modulo) {
        static_assert(std::is_integral<T>::value, "only applicable to integral types");
        return divceil(num, modulo) * modulo;
    }

    uint_t rectmap(uint_t irow, uint_t icol, uint_t ncol);

    void inv_rectmap(uint_t &irow, uint_t &icol, uint_t ncol, uint_t flat);

    uint_t trigmap(uint_t i, uint_t j);

    uint_t npair(uint_t ndim);

    void inv_trigmap(uint_t &i, uint_t &j, uint_t n);

    uint_t strigmap(uint_t i, uint_t j);

    void inv_strigmap(uint_t &i, uint_t &j, uint_t n);

    uint_t nspair(uint_t ndim);

    uint_t factorial(uint_t n);

    uint_t combinatorial(uint_t n, uint_t r);

    uint_t combinatorial_with_repetition(uint_t n, uint_t r);

    template<uint_t n>
    static typename std::enable_if<n == 0ul, uint_t>::type ntup_num(uint_t /*extent*/) {
        return 1ul;
    }

    template<uint_t n>
    static typename std::enable_if<n != 0ul, uint_t>::type ntup_num(uint_t extent) {
        return extent * ntup_num<n - 1>(extent - 1);
    }

    template<uint_t n>
    static uint_t ntup(uint_t extent) {
        return ntup_num<n>(extent) / ntup_num<n>(n);
    }

    template<typename T>
    static std::vector<T> inc(const std::vector<T>& v) {
        auto out = v;
        for (auto& it: out) ++it;
        return out;
    }

    template<typename T>
    static std::vector<T> dec(const std::vector<T>& v) {
        auto out = v;
        for (auto& it: out) --it;
        return out;
    }
}

#endif //M7_UTIL_INTEGER_H
