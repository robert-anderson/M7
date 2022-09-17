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
    std::function<bool(uint_t, uint_t)>
    static make_compare_fn(std::function<typename viewable_t::cview_t(uint_t)> getter_fn, bool asc, bool abs_val) {
        if (asc) {
            if (abs_val)
                return [getter_fn](uint_t i1, uint_t i2) {
                    return std::abs(getter_fn(i1)) < std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](uint_t i1, uint_t i2) {
                    return arith::real(getter_fn(i1)) < arith::real(getter_fn(i2));
                };
        } else {
            if (abs_val)
                return [getter_fn](uint_t i1, uint_t i2) {
                    return std::abs(getter_fn(i1)) > std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](uint_t i1, uint_t i2) {
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
    uintv_t inds(const v_t<T> &v, bool asc, bool abs_val) {
        uintv_t out(v.size());
        std::iota(out.begin(), out.end(), 0);
        if (asc) {
            if (abs_val)
                std::stable_sort(out.begin(), out.end(), [&v](uint_t i, uint_t j) {
                    return std::abs(v[i]) < std::abs(v[j]);
                });
            else
                std::stable_sort(out.begin(), out.end(), [&v](uint_t i, uint_t j) {
                    return arith::real(v[i]) < arith::real(v[j]);
                });
        } else {
            if (abs_val)
                std::stable_sort(out.begin(), out.end(), [&v](uint_t i, uint_t j) {
                    return std::abs(v[i]) > std::abs(v[j]);
                });
            else
                std::stable_sort(out.begin(), out.end(), [&v](uint_t i, uint_t j) {
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
    void inplace(v_t<T> &v, bool asc, bool abs_val) {
        typedef const T &cr_t;
        if (asc) {
            if (abs_val)
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return std::abs(v1) < std::abs(v2);
                });
            else
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return arith::real(v1) < arith::real(v2);
                });
        } else {
            if (abs_val)
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return std::abs(v1) > std::abs(v2);
                });
            else
                std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2) {
                    return arith::real(v1) > arith::real(v2);
                });
        }
    }

    /**
     * apply a naive (copying) reorder operation on a vector determined by the indices in the order vector
     * @tparam T
     *  type of data held by vector being reordered
     * @param in
     *  input vector to be out-of-place reordered
     * @param order
     *  order of the indices of the in vector to appear in the out vector
     * @result
     *  reordered vector
     */
    template<typename T>
    v_t<T> reorder(const v_t<T>& in, const uintv_t& order) {
        v_t<T> out;
        out.reserve(in.size());
        for (auto& i: order) out.push_back(in[i]);
        return out;
    }
}

#endif //M7_SORT_H
