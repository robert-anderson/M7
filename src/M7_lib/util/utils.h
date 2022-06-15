
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

    static constexpr size_t min(size_t i, size_t j) {
        return i < j ? i : j;
    }

    static constexpr size_t max(size_t i, size_t j) {
        return i > j ? i : j;
    }

    /*
     * tests whether all bytes are zero.
     * This is slow, only use for checking that memsets and other zeroing operations
     * have been successful in debug build
     */
    template<typename T>
    static bool is_zero(const T *v) {
        for (size_t ichar = 0ul; ichar < sizeof(T); ++ichar) {
            if (reinterpret_cast<char *>(v)[ichar] != 0) return false;
        }
        return true;
    }

    template<typename T>
    static bool is_zero(const T &v) {
        return is_zero(&v);
    }

    /**
     * catch-all for non-arithmetic types attempts to use output stream << overload. If it doesn't exist, compile time
     * error will be raised
     * @tparam T
     *  type of data to represent as a string
     * @param v
     *  value to represent as string
     * @return
     *  string representation of v
     */
    template<typename T>
    static typename std::enable_if<!std::is_arithmetic<T>::value, std::string>::type
    to_string(const T &v) {
        std::stringstream out;
        out << v;
        return out.str();
    }

    /**
     * catch-all for arithmetic types
     */
    template<typename T>
    static typename std::enable_if<std::is_arithmetic<T>::value, std::string>::type
    to_string(const T &v) {
        if (v == std::numeric_limits<T>::max()) return "inf";
        return std::to_string(v);
    }

    template<typename T>
    static typename std::enable_if<std::is_arithmetic<T>::value, std::string>::type
    to_string(const std::complex<T> &v) {
        if (v == std::numeric_limits<T>::max()) return "inf";
        std::stringstream str;
        str << v;
        return str.str();
    }

    static std::string to_string(const std::vector<std::string> &v) {
        std::string string("[");
        for (const auto &str: v) string += str + " ";
        string += "]";
        return string;
    }

    template<typename T>
    static std::string to_string(const std::vector<T> &v) {
        std::string string("[");
        for (const auto &i: v) string += to_string(i) + " ";
        string += "]";
        return string;
    }

    template<typename T>
    static std::string to_string(const std::vector<std::vector<T>> &v) {
        std::string string("[");
        for (const auto &i: v) string += to_string(i) + " ";
        string += "]";
        return string;
    }

    template<typename T>
    static std::string to_string(const std::stack<T> &v) {
        auto cpy = v;
        std::vector<T> tmp;
        while (!v.empty()) {
            tmp.push_back(cpy.top());
            cpy.pop();
        }
        return to_string(tmp);
    }

    template<typename T>
    static std::string to_string(const std::list<T> &v) {
        auto cpy = v;
        std::vector<T> tmp;
        for (const auto &i: v) tmp.push_back(i);
        return to_string(tmp);
    }

    static std::string to_string(const std::string &str) {
        return "\"" + str + "\"";
    }

    static std::string to_string(bool v) {
        return v ? "true" : "false";
    }

    template<typename T>
    static std::vector<std::string> to_strings(const std::vector<T> &v) {
        std::vector<std::string> out;
        out.reserve(v.size());
        for (auto& i: v) out.push_back(to_string(i));
        return out;
    }

    template<typename T>
    static std::string fp_to_string(const T &v, size_t fp_precision = 6) {
        ASSERT(std::is_floating_point<T>::value);
        std::stringstream tmp;
        tmp << std::scientific << std::setprecision(fp_precision) << v;
        return tmp.str();
    }

    static void pad_string(std::string &str, const size_t num, const char pad_char = ' ') {
        if (num > str.size()) str.insert(0, num - str.size(), pad_char);
    }

    static std::string padded_string(const std::string &str, const size_t num, const char pad_char = ' ') {
        auto tmp = str;
        pad_string(tmp, num, pad_char);
        return tmp;
    }


    template<typename T>
    static std::string num_to_string(const T &entry, size_t padding = 0, size_t fp_precision = 11) {
        std::string result;
        if (std::is_floating_point<T>::value) result = fp_to_string(entry, fp_precision);
        else if (std::is_integral<T>::value) result = std::to_string(entry);
        return result;
    }

    template<typename T>
    static std::string num_to_string(const std::complex<T> &entry, size_t padding = 0, size_t fp_precision = 11) {
        auto tmp_string = fp_to_string(entry.real(), fp_precision) +
                          (entry.imag() < 0 ? "" : "+") + fp_to_string(entry.imag(), fp_precision) + "i";
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }


    template<typename T>
    static void print(typename std::vector<T>::const_iterator begin, typename std::vector<T>::const_iterator end) {
        for (auto iter = begin; iter != end; iter++) {
            std::cout << *iter << " ";
        }
        std::cout << std::endl;
    }

    template<typename T>
    static void print(const std::vector<T> &v) {
        print<T>(v.cbegin(), v.cend());
    }

    template<typename narrow_t, typename wide_t>
    static narrow_t safe_narrow(const wide_t &wide) {
        static_assert(std::is_convertible<wide_t, narrow_t>::value, "incompatible types");
        static_assert(sizeof(wide_t) >= sizeof(narrow_t), "wide type must be at least as long as narrow type");
        ASSERT(static_cast<wide_t>(static_cast<narrow_t>(wide)) == wide); // narrowing loses information
        return static_cast<narrow_t>(wide);
    }

    template<typename narrow_t, typename wide_t>
    static std::vector<narrow_t> safe_narrow(const std::vector<wide_t> &wides) {
        std::vector<narrow_t> narrows;
        narrows.reserve(wides.size());
        for (auto &it : wides) narrows.push_back(utils::safe_narrow<defs::mpi_count>(it));
        return narrows;
    }

    template<typename T1, typename T2>
    static void convert(const std::vector<T1>& v1, std::vector<T2>& v2){
        static_assert(std::is_convertible<T1, T2>::value, "incompatible types");
        v2.clear();
        v2.reserve(v1.size());
        for (auto& i: v1) v2.push_back(i);
    }
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

namespace sort_utils {

    template<typename viewable_t>
    std::function<bool(const size_t &, const size_t &)>
    static make_compare_fn(std::function<typename viewable_t::cview_t(const size_t &)> getter_fn, bool max,
                           bool abs_val) {
        if (max) {
            if (abs_val)
                return [getter_fn](const size_t &i1, const size_t &i2) {
                    return std::abs(getter_fn(i1)) >= std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](const size_t &i1, const size_t &i2) {
                    return getter_fn(i1) >= getter_fn(i2);
                };
        } else {
            if (abs_val)
                return [getter_fn](const size_t &i1, const size_t &i2) {
                    return std::abs(getter_fn(i1)) <= std::abs(getter_fn(i2));
                };
            else
                return [getter_fn](const size_t &i1, const size_t &i2) {
                    return getter_fn(i1) <= getter_fn(i2);
                };
        }
    }

    /**
     * @tparam T
     *  element type
     * @param v
     *  vector to sort
     * @param max
     *  if true, the maximum value indices appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes in v
     * @return
     *  index vector which would sort v to the desired ordering
     */
    template<typename T>
    defs::inds inds(const std::vector<T>& v, bool max, bool abs_val) {
        defs::inds out(v.size());
        std::iota(out.begin(), out.end(), 0);
        if (max) {
            if (abs_val) std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return std::abs(v[i]) < std::abs(v[j]);});
            else std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i] < v[j];});
        } else {
            if (abs_val) std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return std::abs(v[i]) > std::abs(v[j]);});
            else std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i] > v[j];});
        }
        return out;
    }

    /**
     * @tparam T
     *  floating point element type
     * @param v
     *  vector to sort
     * @param max
     *  if true, the maximum value indices appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the complex numbers in v, else compare only the real parts
     * @return
     *  index vector which would sort v to the desired ordering
     */
    template<typename T>
    defs::inds inds(const std::vector<std::complex<T>>& v, bool max, bool abs_val) {
        defs::inds out(v.size());
        std::iota(out.begin(), out.end(), 0);
        if (max) {
            if (abs_val) std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return std::abs(v[i]) < std::abs(v[j]);});
            else std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i].real() < v[j].real();});
        } else {
            if (abs_val) std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){
                    return std::abs(v[i]) > std::abs(v[j]);});
            else std::sort(out.begin(), out.end(), [&v](size_t i, size_t j){return v[i].real() > v[j].real();});
        }
        return out;
    }

    /**
     * @tparam T
     *  floating point type
     * @param v
     *  vector to sort
     * @param max
     *  if true, the maximum values appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes in v
     */
    template<typename T>
    void inplace(std::vector<T>& v, bool max, bool abs_val) {
        typedef const T& cr_t;
        if (max) {
            if (abs_val) std::sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return std::abs(v1) < std::abs(v2);});
            else std::sort(v.begin(), v.end(), [&v](cr_t v1, cr_t v2){ return v1 < v2;});
        } else {
            if (abs_val) std::sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return std::abs(v1) > std::abs(v2);});
            else std::sort(v.begin(), v.end(), [&v](cr_t v1, cr_t v2){ return v1 > v2;});
        }
    }

    /**
     * @tparam T
     *  component floating point type
     * @param v
     *  vector to sort
     * @param max
     *  if true, the maximum values appear first in the output
     * @param abs_val
     *  if true, sort according to the magnitudes of the complex numbers in v, else compare only the real parts
     */
    template<typename T>
    void inplace(std::vector<std::complex<T>>& v, bool max, bool abs_val) {
        typedef const std::complex<T>& cr_t;
        if (max) {
            if (abs_val) std::sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return std::abs(v1) > std::abs(v2);});
            else std::sort(v.begin(), v.end(), [&v](cr_t v1, cr_t v2){ return v1.real() > v2.real();});
        } else {
            if (abs_val) std::sort(v.begin(), v.end(), [](cr_t v1, cr_t v2){ return std::abs(v1) < std::abs(v2);});
            else std::sort(v.begin(), v.end(), [&v](cr_t v1, cr_t v2){ return v1.real() < v2.real();});
        }
    }

    template<typename T>
    void inplace(std::vector<std::complex<T>>& v, bool max) {
        inplace(v, max, true);
    }

}


namespace tags {
    template<typename T>
    struct Type {};

    template<size_t n>
    struct Int {
        static constexpr size_t value(){return n;}
    };
}

template<typename T>
static std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << utils::to_string(v);
    return os;
}


#endif //M7_UTILS_H
