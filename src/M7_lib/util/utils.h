
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_UTILS_H
#define M7_UTILS_H

#include <M7_lib/defs.h>
#include <vector>
#include <array>
#include <stack>
#include <complex>
#include <iostream>
#include <cmath>
#include <limits>
#include <numeric>
#include <x86intrin.h>
#include <climits>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <memory>
#include <functional>
#include <list>

namespace utils {
}

namespace prob_utils {

    static void normalize(std::vector<defs::prob_t> &v, defs::prob_t norm = 1.0) {
        auto tot = std::accumulate(v.begin(), v.end(), 0.0);
        auto fac = norm / tot;
        for (auto &i:v) i *= fac;
    }

    static void rectify(std::vector<defs::prob_t> &v, defs::prob_t min) {
        for (auto &prob: v) if (prob < min) prob = min;
        normalize(v);
    }

    static defs::prob_t linear_bias_prob(const size_t &n, const size_t &i) {
        return (2 * i + 1) / defs::prob_t(n * n);
    }
}

namespace stat_utils {

    template<typename T>
    std::pair<T, T> mean_std(
            typename std::vector<T>::const_iterator begin,
            typename std::vector<T>::const_iterator end) {
        T mean = 0.0;
        T sq_mean = 0.0;
        for (auto i = begin; i != end; i++) {
            mean += *i;
            sq_mean += (*i) * (*i);
        }
        mean /= std::distance(begin, end);
        sq_mean /= std::distance(begin, end);
        return {mean, std::sqrt(std::abs(sq_mean - mean * mean))};
    }

    template<typename T>
    std::pair<T, T> product(const std::pair<T, T> &a, const std::pair<T, T> &b) {
        /*
         * combine statistics in quadrature for product of random variables
         * a and b assuming no covariance
         */
        return {
                a.first * b.first, std::abs(a.first * b.first) *
                                   std::sqrt(std::pow(a.second / a.first, 2.0) + std::pow(b.second / b.first, 2.0))
        };
    }

    template<typename T>
    std::pair<T, T> quotient(const std::pair<T, T> &a, const std::pair<T, T> &b) {
        /*
         * combine statistics in quadrature for quotient of random variables
         * a and b assuming no covariance
         */
        return {
                a.first / b.first, std::abs(a.first / b.first) *
                                   std::sqrt(std::pow(a.second / a.first, 2.0) + std::pow(b.second / b.first, 2.0))
        };
    }

}

namespace float_utils {
    template<typename T>
    bool is_integral(const T &v) {
        static_assert(std::is_floating_point<T>::value, "T must be floating point");
        return ceil(v) == v;
    }
}

namespace complex_utils {
    template<typename T, typename U>
    void set_imag_part(T &v, const U &imag) {}

    template<typename T, typename U>
    void set_imag_part(std::complex<T> &v, const U &imag) { v.imag(T(imag)); }

    template<typename T>
    static std::complex<double> normal_from_polar(const T &arg) {
        /*
         * arg is in radian
         */
        return {std::cos(arg), std::sin(arg)};
    }

    template<typename T>
    static std::complex<double> normal_from_xy(const T &x, const T &y) {
        return std::complex<double>(x, y) / std::sqrt(x * x + y * y);
    }

    static std::complex<double> normal_from_sector(const size_t &isector, const size_t &nsector) {
        return normal_from_polar(consts::two_pi * double(isector) / double(nsector));
    }

    template<typename T>
    static void combine(const std::vector<T>& real, const std::vector<T>& imag, std::vector<std::complex<T>>& v){
        auto n = std::min(real.size(), imag.size());
        v.clear();
        v.reserve(n);
        for (size_t i=0ul; i<n; ++i) v.push_back({real[i], imag[i]});
    }
}

namespace mem_utils {

    template<typename T, typename... Args>
    static std::unique_ptr<T> make_unique(Args &&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    static void print_cmp(char *c1, char *c2, size_t n) {
        for (size_t i = 0; i < n; ++i) {
            std::cout << int(c1[i]) << " " << int(c2[i]) << std::endl;
        }
        std::cout << std::endl;
    }
}

#endif //M7_UTILS_H
