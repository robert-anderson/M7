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

    /**
     * @tparam T
     *  floating point element type
     * @param order
     *  output indices which would put the data in the required order
     * @param v
     *  data to sort
     * @param number of elements in data
     *  data to sort
     * @param asc
     *  if true, the lowest value appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the numbers in v, else compare only the real parts
     * @return
     *  index vector which would sort v to the desired ordering
     */
    template<typename T>
    void inds(uintv_t& order, const T* v, uint_t nelement, bool asc, bool abs_val) {
        order.resize(nelement);
        std::iota(order.begin(), order.end(), 0);
        if (asc) {
            if (abs_val)
                std::stable_sort(order.begin(), order.end(), [&v](uint_t i, uint_t j) {
                    return std::abs(v[i]) < std::abs(v[j]);
                });
            else
                std::stable_sort(order.begin(), order.end(), [&v](uint_t i, uint_t j) {
                    return arith::real(v[i]) < arith::real(v[j]);
                });
        } else {
            if (abs_val)
                std::stable_sort(order.begin(), order.end(), [&v](uint_t i, uint_t j) {
                    return std::abs(v[i]) > std::abs(v[j]);
                });
            else
                std::stable_sort(order.begin(), order.end(), [&v](uint_t i, uint_t j) {
                    return arith::real(v[i]) > arith::real(v[j]);
                });
        }
    }


    template<typename T>
    uintv_t inds(const T* v, uint_t nelement, bool asc, bool abs_val) {
        uintv_t out;
        inds(out, v, nelement, asc, abs_val);
        return out;
    }


    template<typename T>
    uintv_t inds(const v_t<T> &v, bool asc, bool abs_val) {
        return inds(v.data(), v.size(), asc, abs_val);
    }

    /**
     * @tparam T
     *  component floating point type
     * @param v
     *  data to sort
     * @param nelement
     *  number of elements in data
     * @param asc
     *  if true, the lowest value appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the numbers in v, else compare only the real parts
     */
    template<typename T>
    void inplace(T* v, uint_t nelement, bool asc, bool abs_val) {
        typedef const T &cr_t;
        if (asc) {
            if (abs_val)
                std::stable_sort(v, v + nelement, [](cr_t v1, cr_t v2) {
                    return std::abs(v1) < std::abs(v2);
                });
            else
                std::stable_sort(v, v + nelement, [](cr_t v1, cr_t v2) {
                    return arith::real(v1) < arith::real(v2);
                });
        } else {
            if (abs_val)
                std::stable_sort(v, v + nelement, [](cr_t v1, cr_t v2) {
                    return std::abs(v1) > std::abs(v2);
                });
            else
                std::stable_sort(v, v + nelement, [](cr_t v1, cr_t v2) {
                    return arith::real(v1) > arith::real(v2);
                });
        }
    }


    template<typename T>
    void inplace(v_t<T> &v, bool asc, bool abs_val) {
        inplace(v.data(), v.size(), asc, abs_val);
    }

    /**
     * apply a naive (copying) reorder operation on a vector determined by the indices in the order vector
     * @param data
     *  data to be in-place sorted
     * @param element_size
     *  size of each element in data
     * @param order
     *  element indices of initial data in the final data buffer
     */
    void reorder(void* data, uint_t element_size, const uintv_t& order);

    template<typename T>
    void reorder(T* data, const uintv_t& order) {
        reorder(data, sizeof(T), order);
    }

    template<typename T>
    v_t<T> reorder(const v_t<T>& in, const uintv_t& order) {
        v_t<T> out(order.size());
        reorder(out.data(), order);
        return out;
    }
}

#endif //M7_SORT_H
