//
// Created by rja on 13/06/22.
//

#include "Integer.h"
#include "M7_lib/parallel/MPIAssert.h"

uint_t integer::rectmap(uint_t irow, uint_t icol, uint_t ncol) {
    // rectangular map
    /*
     * n         i,j
     * -----------------
     * 0 1 2 3   0,0 0,1 0,2 0,3
     * 4 5 6 7   1,0 1,1 1,2 1,3
     */
    return irow * ncol + icol;
}

void integer::inv_rectmap(uint_t &irow, uint_t &icol, uint_t ncol, uint_t flat) {
    irow = flat / ncol;
    icol = flat - irow * ncol;
}

uint_t integer::trigmap(uint_t i, uint_t j) {
    ASSERT(i >= j);
    /*
     * n        i,j
     * -----------------
     * 0        0,0
     * 1 2      1,0 1,1
     * 3 4 5    2,0 2,1 2,2
     * 6 7 8 9  3,0 3,1 3,2 3,3
     */
    return (i * (i + 1)) / 2 + j;
}

uint_t integer::npair(uint_t ndim) {
    return trigmap(ndim, 0);
}

void integer::inv_trigmap(uint_t &i, uint_t &j, uint_t n) {
    i = (uint_t) ((std::sqrt(1 + 8 * (double) n) - 1) / 2);
    j = n - (i * (i + 1)) / 2;
}

uint_t integer::strigmap(uint_t i, uint_t j) {
    ASSERT(i > j);
    // strict triangular map i>j
    /*
     * n        i,j
     * -----------------
     * 0        1,0
     * 1 2      2,0 2,1
     * 3 4 5    3,0 3,1 3,2
     * 6 7 8 9  4,0 4,1 4,2 4,3
     */
    return (i * (i - 1)) / 2 + j;
}

void integer::inv_strigmap(uint_t &i, uint_t &j, uint_t n) {
    i = (uint_t) ((std::sqrt(1 + 8 * (double) n) + 1) / 2);
    j = n - (i * (i - 1)) / 2;
}

uint_t integer::nspair(uint_t ndim) {
    if (!ndim) return 0ul;
    return strigmap(ndim, 0);
}

uint_t integer::factorial(uint_t n) {
    DEBUG_ASSERT_LT(n, ((uint_t) -1) / 2, "invalid n");
    uint_t out = 1ul;
    if (n < 1) return 1ul;
    for (uint_t i = 1ul; i <= n; ++i) out *= i;
    return out;
}

uint_t integer::combinatorial(uint_t n, uint_t r) {
    /*
     * n choose r = n! / ((n-r)!r!)
     * compute numerator and denominator simultaneously whenever an
     * exact quotient can be computed to avoid premature overflow
     */
    if (r > n) return 0ul;
    if (r == 0) return 1ul;
    if (n == 1) return 1ul;
    if (r == n) return 1ul;

    uint_t out = 1ul;
    uint_t ni = 0ul;
    uint_t ri = 0ul;
    while (1) {
        if (ri < r && out % (r - ri) == 0) {
            out /= r - (ri++);
        } else out *= n - (ni++);
        DEBUG_ASSERT_LE(ni, r, "overflow occurred in exact combinatorial");
        if (ri == r && ni == r) return out;
    }
}

uint_t integer::combinatorial_with_repetition(uint_t n, uint_t r) {
    return combinatorial(n + r - 1, r);
}

uint_t integer::gcd(uint_t a, uint_t b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}

uint_t integer::lcm(uint_t a, uint_t b) {
    return (a*b)/gcd(a, b);
}

uint_t integer::lcm_le(uint_t n) {
    if (!n) return 0;
    uint_t out = 1;
    for (uint_t i = 2; i <= n; i++) out = (i * out) / gcd(i, out);
    return out;
}
