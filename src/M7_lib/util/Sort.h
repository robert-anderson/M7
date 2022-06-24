//
// Created by rja on 15/06/22.
//

#ifndef M7_SORT_H
#define M7_SORT_H

#include <algorithm>
#include <numeric>
#include "M7_lib/defs.h"
#include "Arith.h"

namespace sort {

    template<typename viewable_t>
    std::function<bool(size_t, size_t)>
    static make_compare_fn(std::function<typename viewable_t::cview_t(size_t)> getter_fn, bool asc, bool abs_val) {
        if (asc) {
            if (abs_val)
                return [getter_fn](size_t i1, size_t i2) {
                    return std::abs(getter_fn(i1)) < std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](size_t i1, size_t i2) {
                    return arith::real(getter_fn(i1)) < arith::real(getter_fn(i2));
                };
        } else {
            if (abs_val)
                return [getter_fn](size_t i1, size_t i2) {
                    return std::abs(getter_fn(i1)) > std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](size_t i1, size_t i2) {
                    return arith::real(getter_fn(i1)) > arith::real(getter_fn(i2));
                };
        }
    }

    /**
     * @tparam T
     *  floating point element type
     * @param v
     *  vector to sort
     * @param asc
     *  if true, the lowest value appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the numbers in v, else compare only the real parts
     * @return
     *  index vector which would sort v to the desired ordering
     */
    template<typename T>
    defs::uintv_t inds(const std::vector<T> &v, bool asc, bool abs_val) {
        defs::uintv_t out(v.size());
        std::iota(out.begin(), out.end(), 0);
        if (asc) {
            if (abs_val)
                std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j) {
                    return std::abs(v[i]) < std::abs(v[j]);
                });
            else
                std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j) {
                    return arith::real(v[i]) < arith::real(v[j]);
                });
        } else {
            if (abs_val)
                std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j) {
                    return std::abs(v[i]) > std::abs(v[j]);
                });
            else
                std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j) {
                    return arith::real(v[i]) > arith::real(v[j]);
                });
        }
        return out;
    }

    /**
     * @tparam T
     *  component floating point type
     * @param v
     *  vector to sort
     * @param asc
     *  if true, the lowest value appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the numbers in v, else compare only the real parts
     */
    template<typename T>
    void inplace(std::vector<T> &v, bool asc, bool abs_val) {
        typedef const T &cr_t;
        if (asc) {
            if (abs_val)
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return std::abs(v1) > std::abs(v2);
                });
            else
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return arith::real(v1) > arith::real(v2);
                });
        } else {
            if (abs_val)
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return std::abs(v1) < std::abs(v2);
                });
            else
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return arith::real(v1) < arith::real(v2);
                });
        }
    }
}

#endif //M7_SORT_H
