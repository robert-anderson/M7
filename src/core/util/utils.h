//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_UTILS_H
#define M7_UTILS_H

#include "defs.h"
#include <vector>
#include <complex>
#include <iostream>
#include <cmath>
#include <limits>
#include <numeric>
#include <x86intrin.h>
#include <climits>
#include <iomanip>

namespace utils {
    template<typename T>
    std::string to_string(T v, size_t n) {
        std::string string("[");
        for (size_t i = 0ul; i < n; ++i) {
            string += std::to_string(v[i]) + " ";
        }
        string += "]";
        return string;
    }

    template<typename T>
    std::string to_string(T v) {
        return to_string(v, v.size());
    }

    template<typename T>
    std::string fp_to_string(const T &v, size_t fp_precision = 6) {
        ASSERT(std::is_floating_point<T>::value);
        std::stringstream tmp;
        tmp << std::scientific << std::setprecision(fp_precision) << v;
        return tmp.str();
    }

    template<typename T>
    std::string num_to_string(const T &entry, size_t padding = 0, size_t fp_precision = 6) {
        std::string result;
        if (std::is_floating_point<T>::value) result = fp_to_string(entry, fp_precision);
        else if (std::is_integral<T>::value) result = std::to_string(entry);
        auto decimal_length = std::numeric_limits<T>::digits10;
        //assert(result.size()<=decimal_length);
        //result.insert(result.begin(), padding + decimal_length - result.size(), ' ');
        return result;
    }

    template<typename T>
    std::string num_to_string(const std::complex<T> &entry, size_t padding = 0, size_t fp_precision = 6) {
        auto tmp_string = fp_to_string(entry.real(), fp_precision) +
                          (entry.imag() < 0 ? "" : "+") + fp_to_string(entry.imag(), fp_precision) + "i";
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }


    template<typename T>
    void print(typename std::vector<T>::const_iterator begin, typename std::vector<T>::const_iterator end) {
        for (auto iter = begin; iter != end; iter++) {
            std::cout << *iter << " ";
        }
        std::cout << std::endl;
    }

    template<typename T>
    void print(const std::vector<T> &v) {
        print<T>(v.cbegin(), v.cend());
    }

}

namespace integer_utils {

    template<typename T>
    static typename std::enable_if<std::is_integral<T>::value, T>::type
    divceil(const T &num, const T &denom) {
        return num % denom ? num / denom + 1 : num / denom;
    }

    size_t rectmap(const size_t &irow, const size_t &icol, const size_t &ncol);

    void inv_rectmap(size_t &irow, size_t &icol, const size_t &ncol, const size_t &flat);

    size_t trigmap(const size_t &i, const size_t &j);

    size_t npair(const size_t &ndim);

    void inv_trigmap(size_t &i, size_t &j, const size_t &n);

    size_t strigmap(const size_t &i, const size_t &j);

    void inv_strigmap(size_t &i, size_t &j, const size_t &n);

    size_t nspair(const size_t &ndim);

    size_t factorial(const size_t &n);

    size_t combinatorial(const size_t &n, const size_t &r);
}

namespace bit_utils {
    template<typename T>
    static inline void clr(T &x, size_t i) {
        x &= ~((T) 1ul << i);
    }

    template<typename T>
    static inline void set(T &x, size_t i) {
        x |= ((T) 1ul << i);
    }

    template<typename T>
    static inline bool get(const T &x, size_t i) {
        return (x >> i) & T(1ul);
    }

    template<typename T>
    static inline size_t next_setbit(T &work);

    template<>
    inline size_t next_setbit(unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    template<>
    inline size_t next_setbit(unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    template<>
    inline size_t next_setbit(unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        size_t result = __tzcnt_u32(work);
        bit_utils::clr(work, result);
        return result;
    }

    template<typename T>
    static inline size_t nsetbit(const T &work);

    template<>
    inline size_t nsetbit(const unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    template<>
    inline size_t nsetbit(const unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    template<>
    inline size_t nsetbit(const unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        return _popcnt32(work);
    }

    template<typename T>
    T truncate(T &v, size_t n) {
        const auto nbit = sizeof(T) * CHAR_BIT - n;
        return (v << nbit) >> nbit;
    }

}


namespace string_utils {
    static std::string join(const std::vector<std::string> &words, const std::string &divider, const bool &bookends) {
        std::string out{""};
        if (bookends) out += divider;
        for (size_t i = 0ul; i < words.size() - 1; ++i) {
            out += words[i] + divider;
        }
        out += words[words.size() - 1];
        if (bookends) out += divider;
        return out;
    }

    static std::string join(const std::vector<std::string> &words, const std::string &divider) {
        return join(words, divider, false);
    }

    static std::string join(const std::vector<std::string> &words, const bool &bookends) {
        return join(words, " ", bookends);
    }

    static std::string join(const std::vector<std::string> &words) {
        return join(words, " ", false);
    }

    static std::string
    join(const std::string &word, const size_t &nrepeat, const std::string &divider, const bool &bookends) {
        return join(std::vector<std::string>(nrepeat, word), divider, bookends);
    }

    static std::string join(const std::string &word, const size_t &nrepeat, const std::string &divider) {
        return join(word, nrepeat, divider, false);
    }

    static std::string join(const std::string &word, const size_t &nrepeat) {
        return join(word, nrepeat, " ", false);
    }
}

namespace prob_utils {
    template<typename T>
    void normalize(std::vector<T> &v, T norm = T(1)) {
        T tot = std::accumulate(v.begin(), v.end(), T(0));
        T fac = norm / tot;
        for (auto &i:v) i *= fac;
    }
}

namespace stat_utils {

    template<typename T>
    std::pair<T, T> mean_std(
            typename std::vector<T>::const_iterator begin,
            typename std::vector<T>::const_iterator end){
        T mean = 0.0;
        T sq_mean = 0.0;
        for (auto i = begin; i != end; i++) {
            mean += *i;
            sq_mean += (*i)*(*i);
        }
        mean /= std::distance(begin, end);
        sq_mean /= std::distance(begin, end);
        return {mean, std::sqrt(std::abs(sq_mean-mean*mean))};
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

#endif //M7_UTILS_H
