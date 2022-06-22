//
// Created by rja on 13/06/22.
//

#include "Integer.h"


size_t utils::integer::rectmap(size_t irow, size_t icol, size_t ncol) {
    // rectangular map
    /*
     * n         i,j
     * -----------------
     * 0 1 2 3   0,0 0,1 0,2 0,3
     * 4 5 6 7   1,0 1,1 1,2 1,3
     */
    return irow * ncol + icol;
}

void utils::integer::inv_rectmap(size_t &irow, size_t &icol, size_t ncol, size_t flat) {
    irow = flat / ncol;
    icol = flat - irow * ncol;
}

size_t utils::integer::trigmap(size_t i, size_t j) {
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

size_t utils::integer::npair(size_t ndim) {
    return trigmap(ndim, 0);
}

void utils::integer::inv_trigmap(size_t &i, size_t &j, size_t n) {
    i = (size_t) ((std::sqrt(1 + 8 * (double) n) - 1) / 2);
    j = n - (i * (i + 1)) / 2;
}

size_t utils::integer::strigmap(size_t i, size_t j) {
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

void utils::integer::inv_strigmap(size_t &i, size_t &j, size_t n) {
    i = (size_t) ((std::sqrt(1 + 8 * (double) n) + 1) / 2);
    j = n - (i * (i - 1)) / 2;
}

size_t utils::integer::nspair(size_t ndim) {
    if (!ndim) return 0ul;
    return strigmap(ndim, 0);
}

size_t utils::integer::factorial(size_t n) {
    ASSERT(n < ((size_t) -1) / 2);
    size_t out = 1ul;
    if (n < 1) return 1ul;
    for (size_t i = 1ul; i <= n; ++i) out *= i;
    return out;
}

size_t utils::integer::combinatorial(size_t n, size_t r) {
    /*
     * n choose r = n! / ((n-r)!r!)
     * compute numerator and denominator simultaneously whenever an
     * exact quotient can be computed to avoid premature overflow
     */
    ASSERT(n >= r);
    if (r == 0) return 1ul;
    if (n == 1) return 1ul;
    if (r == n) return 1ul;

    size_t out = 1ul;
    size_t ni = 0ul;
    size_t ri = 0ul;
    while (1) {
        if (ri < r && out % (r - ri) == 0) {
            out /= r - (ri++);
        } else out *= n - (ni++);
        ASSERT(ni <= r); // overflow occurred.
        if (ri == r && ni == r) return out;
    }
}

size_t utils::integer::combinatorial_with_repetition(size_t n, size_t r) {
    return combinatorial(n + r - 1, r);
}