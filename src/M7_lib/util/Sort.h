//
// Created by rja on 15/06/22.
//

#ifndef M7_SORT_H
#define M7_SORT_H

#include "utils.h"

namespace utils {

    namespace sort {

        template<typename viewable_t>
        std::function<bool(size_t, size_t)>
        static make_compare_fn(std::function<typename viewable_t::cview_t(size_t )> getter_fn, bool asc, bool abs_val) {
            if (asc) {
                if (abs_val)
                    return [getter_fn](size_t i1, size_t i2) {
                        return std::abs(getter_fn(i1)) < std::abs(getter_fn(i2));
                    };
                else
                    return [getter_fn](size_t i1, size_t i2) {
                        return consts::real(getter_fn(i1)) < consts::real(getter_fn(i2));
                    };
            } else {
                if (abs_val)
                    return [getter_fn](size_t i1, size_t i2) {
                        return std::abs(getter_fn(i1)) > std::abs(getter_fn(i2));
                    };
                else
                    return [getter_fn](size_t i1, size_t i2) {
                        return consts::real(getter_fn(i1)) > consts::real(getter_fn(i2));
                    };
            }
        }

        /**
         * @tparam T
         *  element type
         * @param v
         *  vector to sort
         * @param asc
         *  if true, the lowest value appear first in the output
         * @param abs_val
         *  if true, sort according to the magnitudes in v
         * @return
         *  index vector which would sort v to the desired ordering
         */
        template<typename T>
        defs::inds inds(const std::vector<T>& v, bool asc, bool abs_val) {
            defs::inds out(v.size());
            std::iota(out.begin(), out.end(), 0);
            if (asc) {
                if (abs_val) std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                        return std::abs(v[i]) < std::abs(v[j]);});
                else std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i] < v[j];});
            } else {
                if (abs_val) std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                        return std::abs(v[i]) > std::abs(v[j]);});
                else std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i] > v[j];});
            }
            return out;
        }

        /**
         * @tparam T
         *  floating point element type
         * @param v
         *  vector to sort
         * @param asc
         *  if true, the lowest value appear first in the output
         * @param abs_val
         *  if true, sort according to the magnitudes of the complex numbers in v, else compare only the real parts
         * @return
         *  index vector which would sort v to the desired ordering
         */
        template<typename T>
        defs::inds inds(const std::vector<std::complex<T>>& v, bool asc, bool abs_val) {
            defs::inds out(v.size());
            std::iota(out.begin(), out.end(), 0);
            if (asc) {
                if (abs_val) std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                        return std::abs(v[i]) < std::abs(v[j]);});
                else std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return v[i].real() < v[j].real();});
            } else {
                if (abs_val) std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                        return std::abs(v[i]) > std::abs(v[j]);});
                else std::stable_sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return v[i].real() > v[j].real();});
            }
            return out;
        }

        /**
         * @tparam T
         *  floating point type
         * @param v
         *  vector to sort
         * @param asc
         *  if true, the lowest value appear first in the output
         * @param abs_val
         *  if true, sort according to the magnitudes in v
         */
        template<typename T>
        void inplace(std::vector<T>& v, bool asc, bool abs_val) {
            typedef const T& cr_t;
            if (asc) {
                if (abs_val) std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return std::abs(v1) < std::abs(v2);});
                else std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return v1 < v2;});
            } else {
                if (abs_val) std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return std::abs(v1) > std::abs(v2);});
                else std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return v1 > v2;});
            }
        }

        /**
         * @tparam T
         *  component floating point type
         * @param v
         *  vector to sort
         * @param asc
         *  if true, the lowest value appear first in the output
         * @param abs_val
         *  if true, sort according to the magnitudes of the complex numbers in v, else compare only the real parts
         */
        template<typename T>
        void inplace(std::vector<std::complex<T>>& v, bool asc, bool abs_val) {
            typedef const std::complex<T>& cr_t;
            if (asc) {
                if (abs_val) std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return std::abs(v1) > std::abs(v2);});
                else std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return v1.real() > v2.real();});
            } else {
                if (abs_val) std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return std::abs(v1) < std::abs(v2);});
                else std::stable_sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){
                    return v1.real() < v2.real();});
            }
        }
    }
}

#endif //M7_SORT_H
