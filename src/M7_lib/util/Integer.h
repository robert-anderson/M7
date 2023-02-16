//
// Created by rja on 13/06/22.
//

#ifndef M7_UTIL_INTEGER_H
#define M7_UTIL_INTEGER_H

#include "M7_lib/defs.h"
#include "M7_lib/parallel/MPIAssert.h"

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
     * @param nitem
     *  total number of items to be distributed among nbin bins
     * @param ibin
     *  index of the specific bin
     * @param nbin
     *  number of bins
     * @return
     *  number of items in the ibin-indexed bin
     */
    static constexpr uint_t evenly_shared_count(uint_t nitem, uint_t ibin, uint_t nbin) {
        return nitem / nbin + (ibin < (nitem % nbin));
    }

    /**
     * @param nitem
     *  total number of items to be distributed among nbin bins
     * @param ibin
     *  index of the specific bin
     * @param nbin
     *  number of bins
     * @return
     *  number of items in the bins with index less than ibin
     */
    static constexpr uint_t evenly_shared_offset(uint_t nitem, uint_t ibin, uint_t nbin) {
        return min(ibin, nitem % nbin) + ibin * (nitem / nbin);
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

    /**
     * rectangular map
     * n         i,j
     * -----------------
     * 0 1 2 3   0,0 0,1 0,2 0,3
     * 4 5 6 7   1,0 1,1 1,2 1,3
     */
    static uint_t rectmap(uint_t irow, uint_t icol, uint_t ncol) {
        return irow * ncol + icol;
    }

    /**
     * inverse of the rectangular map
     */
    static void inv_rectmap(uint_t &irow, uint_t &icol, uint_t ncol, uint_t flat) {
        irow = flat / ncol;
        icol = flat - irow * ncol;
    }

    /**
     * triangular map
     * n        i,j
     * -----------------
     * 0        0,0
     * 1 2      1,0 1,1
     * 3 4 5    2,0 2,1 2,2
     * 6 7 8 9  3,0 3,1 3,2 3,3
     */
    static uint_t trigmap(uint_t i, uint_t j) {
        DEBUG_ASSERT_GE(i, j, "incorrectly ordered args");
        return (i * (i + 1)) / 2 + j;
    }

    /**
     * reorders args if necessary before delegating trigmap
     */
    static uint_t trigmap_unordered(uint_t i, uint_t j) {
        return i >= j ? trigmap(i, j) : trigmap(j, i);
    }

    uint_t npair(uint_t ndim);

    /**
     * inverse of the triangular map
     */
    static void inv_trigmap(uint_t &i, uint_t &j, uint_t n) {
        i = (uint_t) ((std::sqrt(1 + 8 * (double) n) - 1) / 2);
        j = n - (i * (i + 1)) / 2;
    }

    /**
     * strict triangular map
     * n        i,j
     * -----------------
     * 0        1,0
     * 1 2      2,0 2,1
     * 3 4 5    3,0 3,1 3,2
     * 6 7 8 9  4,0 4,1 4,2 4,3
     */
    static uint_t strigmap(uint_t i, uint_t j) {
        DEBUG_ASSERT_GT(i, j, "incorrectly ordered args");
        return (i * (i - 1)) / 2 + j;
    }

    /**
     * reorders args if necessary before delegating strigmap
     */
    static uint_t strigmap_unordered(uint_t i, uint_t j) {
        return i >= j ? strigmap(i, j) : strigmap(j, i);
    }

    /**
     * inverse of the strict triangular map
     */
    static void inv_strigmap(uint_t &i, uint_t &j, uint_t n) {
        i = (uint_t) ((std::sqrt(1 + 8 * (double) n) + 1) / 2);
        j = n - (i * (i - 1)) / 2;
    }

    uint_t nspair(uint_t ndim);

    uint_t factorial(uint_t n);

    uint_t combinatorial(uint_t n, uint_t r);

    uint_t combinatorial_with_repetition(uint_t n, uint_t r);

    uint_t sqrt(uint_t n);

    /**
     * get greatest common divisor of the two given whole numbers
     */
    uint_t gcd(uint_t a, uint_t b);

    /**
     * get least common multiple of the two given whole numbers
     */
    uint_t lcm(uint_t a, uint_t b);

    /**
     * @param n
     *  maximum number in range
     * @return
     *  least common multiple of all whole numbers less than or equal to n
     */
    uint_t lcm_le(uint_t n);

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
    static void shift(v_t<T>& v, bool increment=true) {
        if (increment) for (auto& it: v) ++it;
        else for (auto& it: v) --it;
    }

    template<typename T>
    static v_t<T> shifted(const v_t<T>& v, bool increment=true) {
        auto out = v;
        shift(out, increment);
        return out;
    }

    template<typename T>
    static str_t to_hex_string(T v) {
        static_assert(std::is_integral<T>::value, "only applicable to integral types");
        std::stringstream stream;
        stream << std::hex << v;
        return stream.str();
    }
}

#endif //M7_UTIL_INTEGER_H
