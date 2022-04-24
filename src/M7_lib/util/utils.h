
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

    template<size_t exp, typename T=void>
    static typename std::enable_if<exp == 0, T>::type
    pow(const T &v) {
        return 1ul;
    }

    template<size_t exp, typename T=void>
    static typename std::enable_if<exp != 0, T>::type
    pow(const T &v) {
        return v * pow<exp - 1, T>(v);
    }

    template<size_t n>
    static typename std::enable_if<n == 0ul, size_t>::type
    ntup_num(size_t extent) {
        return 1ul;
    }

    template<size_t n>
    static typename std::enable_if<n != 0ul, size_t>::type
    ntup_num(size_t extent) {
        return extent * ntup_num<n - 1>(extent - 1);
    }

    template<size_t n>
    static size_t ntup(size_t extent) {
        return ntup_num<n>(extent) / ntup_num<n>(n);
    }

    template<typename T1, typename T2>
    static void convert(const std::vector<T1>& v1, std::vector<T2>& v2){
        static_assert(std::is_convertible<T1, T2>::value, "incompatible types");
        v2.clear();
        v2.reserve(v1.size());
        for (auto& i: v1) v2.push_back(i);
    }
}

namespace integer_utils {

    template<typename T>
    static typename std::enable_if<std::is_integral<T>::value, T>::type
    divceil(const T &num, const T &denom) {
        return num % denom ? num / denom + 1 : num / denom;
    }

    template<typename T>
    static typename std::enable_if<std::is_integral<T>::value, T>::type
    round_up(const T &num, const T &modulo) {
        return divceil(num, modulo) * modulo;
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

    size_t combinatorial_with_repetition(const size_t &n, const size_t &r);

    /**
     * exact integer power by recursive squaring
     * @tparam T
     *  integer type
     * @param x
     *  number being exponentiated
     * @param y
     *  exponent
     * @return
     *  x to the power y
     */
    template<typename T>
    T pow(T x, T y){
        static_assert(std::is_integral<T>::value, "template arg must be integral type");
        if (!y) return 1;
        if (y == 1) return x;

        auto tmp = pow(x, y / 2);
        if (y & T(1)) return x * tmp * tmp;
        return tmp * tmp;
    }
}

namespace bit_utils {

    /**
     * masks whose ith elements retain i bits when bitwise and-ed with another value of the same size (4 or 8 bytes)
     */
    static constexpr std::array<uint8_t, 9> c_trunc_mask_8 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff};
    static constexpr std::array<uint16_t, 17> c_trunc_mask_16 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff};
    static constexpr std::array<uint32_t, 33> c_trunc_mask_32 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff,
             0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff,
             0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff, 0x1fffffff,
             0x3fffffff, 0x7fffffff, 0xffffffff};
    static constexpr std::array<uint64_t, 65> c_trunc_mask_64 =
            {0x0, 0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff,
             0x7ff, 0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff,
             0x3ffff, 0x7ffff, 0xfffff, 0x1fffff, 0x3fffff, 0x7fffff,
             0xffffff, 0x1ffffff, 0x3ffffff, 0x7ffffff, 0xfffffff,
             0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff, 0x1ffffffff,
             0x3ffffffff, 0x7ffffffff, 0xfffffffff, 0x1fffffffff,
             0x3fffffffff, 0x7fffffffff, 0xffffffffff, 0x1ffffffffff,
             0x3ffffffffff, 0x7ffffffffff, 0xfffffffffff, 0x1fffffffffff,
             0x3fffffffffff, 0x7fffffffffff, 0xffffffffffff,
             0x1ffffffffffff, 0x3ffffffffffff, 0x7ffffffffffff,
             0xfffffffffffff, 0x1fffffffffffff, 0x3fffffffffffff,
             0x7fffffffffffff, 0xffffffffffffff, 0x1ffffffffffffff,
             0x3ffffffffffffff, 0x7ffffffffffffff, 0xfffffffffffffff,
             0x1fffffffffffffff, 0x3fffffffffffffff, 0x7fffffffffffffff,
             0xffffffffffffffff};

    template<typename T>
    static inline void clr(T &x, const size_t &i) {
        x &= ~(T(1) << T(i));
    }

    template<typename T>
    static inline void set(T &x, const size_t &i) {
        x |= (T(1) << T(i));
    }

    template<typename T>
    static inline bool get(const T &x, const size_t &i) {
        return x & (T(1) << T(i));
    }

    inline size_t next_setbit(unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    inline size_t next_setbit(unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    inline size_t next_setbit(unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        size_t result = __tzcnt_u32(work);
        bit_utils::clr(work, result);
        return result;
    }

    inline size_t nsetbit(const unsigned long long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    inline size_t nsetbit(const unsigned long &work) {
        static_assert(sizeof(work) == 8, "Data length not supported");
        return _popcnt64(work);
    }

    inline size_t nsetbit(const unsigned &work) {
        static_assert(sizeof(work) == 4, "Data length not supported");
        return _popcnt32(work);
    }

    static uint8_t truncate(const uint8_t &v, const size_t& n) {
        return v & c_trunc_mask_8[n];
    }
    static uint16_t truncate(const uint16_t &v, const size_t& n) {
        return v & c_trunc_mask_16[n];
    }
    static uint32_t truncate(const uint32_t &v, const size_t& n) {
        return v & c_trunc_mask_32[n];
    }
    static uint64_t truncate(const uint64_t &v, const size_t& n) {
        return v & c_trunc_mask_64[n];
    }

    static int8_t truncate(const int8_t &v, const size_t& n) {
        return truncate(uint8_t(v), n);
    }
    static int16_t truncate(const int16_t &v, const size_t& n) {
        return truncate(uint16_t(v), n);
    }
    static int32_t truncate(const int32_t &v, const size_t& n) {
        return truncate(uint32_t(v), n);
    }
    static int64_t truncate(const int64_t &v, const size_t& n) {
        return truncate(uint64_t(v), n);
    }

    static void make_range_mask(uint8_t& v, const size_t& ibegin, const size_t& iend){
        v = c_trunc_mask_8[iend] & ~c_trunc_mask_8[ibegin];
    }
    static void make_range_mask(uint16_t& v, const size_t& ibegin, const size_t& iend){
        v = c_trunc_mask_16[iend] & ~c_trunc_mask_16[ibegin];
    }
    static void make_range_mask(uint32_t& v, const size_t& ibegin, const size_t& iend){
        v = c_trunc_mask_32[iend] & ~c_trunc_mask_32[ibegin];
    }
    static void make_range_mask(uint64_t& v, const size_t& ibegin, const size_t& iend){
        v = c_trunc_mask_64[iend] & ~c_trunc_mask_64[ibegin];
    }

    template<typename T>
    static T make_range_mask(const size_t& ibegin, const size_t& iend) {
        T v;
        make_range_mask(v, ibegin, iend);
        return v;
    }

    /**
     * apply an arbitrary bit mask to a single buffer word
     * e.g. for a range mask:
     *  F:    10101010101010110
     *  M:    00000111111000000
     *  F|M:  10101010101010110
     *  F&~M: 10101000000010110
     *
     * @tparam T
     *  buffer type of the bit field (integer of any size)
     * @param v
     *  buffer word
     * @param mask
     *  masking word
     * @param set
     *  true if bits that are set in the mask are to be set, else they are cleared
     */
    template<typename T>
    void apply_mask(T& v, const T& mask, bool set){
        if (set) v |= mask;
        else v&= ~mask;
    }
    /**
     * make the range mask word then apply it
     */
    template<typename T>
    void apply_mask(T& v, const size_t& ibegin, const size_t& iend, bool set){
        T mask;
        make_range_mask(mask, ibegin, iend);
        apply_mask(v, mask, set);
    }

    template<typename T>
    static inline size_t nsetbit_before(const T &v, size_t n) {
        return nsetbit(truncate(v, n));
    }

    template<typename T>
    static std::string to_string(const T &v) {
        std::string tmp;
        for (size_t i = 0ul; i < sizeof(T) * CHAR_BIT; ++i) tmp += get(v, i) ? '1' : '0';
        return tmp;
    }

}


namespace string_utils {
    std::string join(const std::vector<std::string> &words, const std::string &divider);

    std::string join(const std::vector<std::string> &words);

    std::string join(const std::string &word, const size_t &nrepeat, const std::string &divider);

    std::string join(const std::string &word, const size_t &nrepeat);

    std::vector<std::string> split(const std::string &line, char delimiter);

    std::vector<std::string> split(const std::string &line, const std::string &delimiters);

    void split(std::string &line, std::vector<std::string>& tokens, const std::string &delimiters);

    std::string yn(bool t);

    std::string YN(bool t);

    std::string memsize(size_t nbyte);

    std::string boxed(std::string s, size_t padding = 4, char c = '#');


    inline bool is_numeric(const char &c);

    inline bool is_partial_standard_float(const char &c);

    inline bool is_partial_scientific(const char &c);

    inline bool is_divider(const char &c);

    double read_double(const char *&ptr);

    size_t read_unsigned(const char *&ptr);

    int64_t read_signed(const char *&ptr);

    size_t parse_decimal_digit(const char *c);

    std::string plural(size_t i, std::string plu_ending="s", std::string sing_ending="");

    std::string plural(std::string base, size_t i, std::string plu_ending="s", std::string sing_ending="");

    std::string prefix(std::string base, std::string prefix="", char delimiter=' ');

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

#if 0
namespace ci_utils {
    /**
     * @param nelec
     *  total number of electrons in determinant
     * @param ms2
     *  total z component of spin 2*(nalpha-nbeta)
     * @return
     *  number of alpha (spin channel 0) electrons
     */
    static size_t nalpha(const FrmHilbertData& hd){
        auto nelec = hd.m_nelec;
        auto ms2 = hd.m_ms2_restrict;
        size_t spin_odd = std::abs(ms2) % 2;
        ASSERT(nelec % 2 == spin_odd);
        size_t nalpha = nelec / 2 + (std::abs(ms2)) / 2 + spin_odd;
        return ms2 >= 0 ? nalpha : nelec - nalpha;
    }

    static size_t nbeta(const FrmHilbertData& hd){
        return hd.m_nelec - nalpha(hd);
    }

    static size_t fermion_dim(size_t nsite, size_t nelec) {
        return integer_utils::combinatorial(2 * nsite, nelec);
    }

    static size_t fermion_dim(size_t nsite, const FrmHilbertData& hd){
        ASSERT(static_cast<size_t>(ms2_restrict) % 2 == nelec % 2)
        auto na = integer_utils::combinatorial(nsite, nalpha(nelec, ms2_restrict));
        auto nb = integer_utils::combinatorial(nsite, nbeta(nelec, ms2_restrict));
        return na*nb;
    }

    /**
     * @param nmode
     *  number of boson modes
     * @param nboson
     *  total number of bosons if number conserving, else the cutoff occupation for each mode
     * @param number_conserve
     *  true if the hamiltonian is boson number-conserving (no ladder terms)
     * @return
     *  dimension of the bosonic sector of the Hilbert space
     */
    static size_t boson_dim(size_t nmode, size_t nboson, bool number_conserve) {
        if (number_conserve) return integer_utils::combinatorial_with_repetition(nmode, nboson);
        else return std::pow(nboson + 1, nmode);
    }
}
#endif

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
}

namespace tuple_utils {
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each(std::tuple<Tp...> &, FuncT &) // Unused arguments are given no names.
    {}

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each(std::tuple<Tp...> &t, FuncT &f) {
        f(std::get<I>(t));
        for_each<I + 1, FuncT, Tp...>(t, f);
    }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each(const std::tuple<Tp...> &, FuncT &) // Unused arguments are given no names.
    {}

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each(const std::tuple<Tp...> &t, FuncT &f) {
        f(std::get<I>(t));
        for_each<I + 1, FuncT, Tp...>(t, f);
    }


    /*
     * modifiable / modifiable
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &, std::tuple<Tp...> &, FuncT &) // Unused arguments are given no names.
    {}

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &t1, std::tuple<Tp...> &t2, FuncT &f) {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }


    /*
     * modifiable / const
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &, const std::tuple<Tp...> &, FuncT &) // Unused arguments are given no names.
    {}

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &t1, const std::tuple<Tp...> &t2, FuncT &f) {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }

    /*
     * const / const
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(const std::tuple<Tp...> &, const std::tuple<Tp...> &, FuncT &) // Unused arguments are given no names.
    {}

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(const std::tuple<Tp...> &t1, const std::tuple<Tp...> &t2, FuncT &f) {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }
}

namespace nd_utils {
    template<typename T>
    T nelement(const std::vector<T> &v) {
        T out = 1;
        for (const auto &i: v) out *= i;
        return out;
    };
}

namespace tags {
    template<typename T>
    struct Type {};

    template<size_t n>
    struct Int {
        static constexpr size_t value(){return n;}
    };
}

namespace array_utils {
    template<typename T, size_t nind>
    static std::array<T, nind> filled(const T &v) {
        std::array<T, nind> tmp;
        tmp.fill(v);
        return tmp;
    }

    template<typename T, size_t nind>
    static std::vector<T> to_vector(const std::array<T, nind> &array) {
        std::vector<T> tmp;
        tmp.assign(array.cbegin(), array.cend());
        return tmp;
    }
}


/**
 * functions related to "excitation signatures" of connections between many-body basis functions
 */
namespace exsig_utils {

    using namespace defs;

    /**
     * compactly expresses an arbitrary SQ operator product as a single integer given some compile-time constant numbers
     * of bits for each element. e.g. if nbit_exsig_nop_frm = 3 and nbit_exsig_nop_bos = 1, then 2x3+2x1 = 8 bits are
     * required to store a connection excitation level as an exsig (excitation signature) with upto 7 fermion creation
     * operators, 7 fermion annihilation operators and 1 each of boson creation and annihilation operators, this limit
     * should be sufficient for all foreseeable applications, but these bit segment lengths are not hardcoded.
     * @param nfrm_cre
     *  number of fermion creation indices in the SQ operator product
     * @param nfrm_ann
     *  number of fermion annihilation indices in the SQ operator product
     * @param nbos_cre
     *  number of boson creation indices in the SQ operator product
     * @param nbos_ann
     *  number of boson annihilation indices in the SQ operator product
     * @return
     *  the excitation signature
     */
    static constexpr size_t encode(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann) {
        return (nfrm_cre > exsig_nop_mask_frm || nfrm_ann > exsig_nop_mask_frm ||
                nbos_cre > exsig_nop_mask_bos || nbos_ann > exsig_nop_mask_bos) ?
               ~0ul : nfrm_cre | (nfrm_ann << nbit_exsig_nop_frm) |
                      (nbos_cre << (2 * nbit_exsig_nop_frm)) |
                      (nbos_ann << (2 * nbit_exsig_nop_frm + nbit_exsig_nop_bos));
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the number of fermion creation indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nfrm_cre(size_t exsig) {
        return exsig_nop_mask_frm & exsig;
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the number of fermion annihilation indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nfrm_ann(size_t exsig) {
        return exsig_nop_mask_frm & (exsig >> nbit_exsig_nop_frm);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the number of boson creation indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nbos_cre(size_t exsig) {
        return exsig_nop_mask_bos & (exsig >> (2 * nbit_exsig_nop_frm));
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the number of boson annihilation indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nbos_ann(size_t exsig) {
        return exsig_nop_mask_bos & (exsig >> (2 * nbit_exsig_nop_frm + nbit_exsig_nop_bos));
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the total number of fermion indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nfrm(size_t exsig) {
        return decode_nfrm_cre(exsig) + decode_nfrm_ann(exsig);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the total number of boson indices in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nbos(size_t exsig) {
        return decode_nbos_cre(exsig) + decode_nbos_ann(exsig);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the total number of operators of any particle type in the SQ operator product encoded within exsig
     */
    static constexpr size_t decode_nop(size_t exsig) {
        return decode_nfrm(exsig) + decode_nbos(exsig);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  true if the exsig has no boson operators
     */
    static constexpr bool is_pure_frm(size_t exsig) {
        return !(decode_nbos_cre(exsig) || decode_nbos_ann(exsig));
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  true if the exsig has no fermion operators
     */
    static constexpr bool is_pure_bos(size_t exsig) {
        return !(decode_nfrm_cre(exsig) || decode_nfrm_ann(exsig));
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  true if the exsig represents a fermion number-conserving operator product
     */
    static constexpr bool conserves_nfrm(size_t exsig) {
        return decode_nfrm_cre(exsig) == decode_nfrm_ann(exsig);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  true if the exsig represents a boson number-conserving operator product
     */
    static constexpr bool conserves_nbos(size_t exsig) {
        return decode_nbos_cre(exsig) == decode_nbos_ann(exsig);
    }

    static constexpr size_t ncontrib_frm(size_t exsig) {
        return utils::min(decode_nfrm_cre(exsig), decode_nfrm_ann(exsig)) + 1;
    }

    static constexpr size_t ncontrib_bos(size_t exsig) {
        return utils::min(decode_nbos_cre(exsig), decode_nbos_ann(exsig)) + 1;
    }

    static constexpr size_t base_exsig(size_t exsig) {
        return encode(
                decode_nfrm_cre(exsig) - (ncontrib_frm(exsig) - 1),
                decode_nfrm_ann(exsig) - (ncontrib_frm(exsig) - 1),
                decode_nbos_cre(exsig) - (ncontrib_bos(exsig) - 1),
                decode_nbos_ann(exsig) - (ncontrib_bos(exsig) - 1)
        );
    }

    static constexpr size_t add_ops(size_t exsig, size_t nfrm, size_t nbos) {
        return encode(decode_nfrm_cre(exsig) + nfrm, decode_nfrm_ann(exsig) + nfrm,
                      decode_nbos_cre(exsig) + nbos, decode_nbos_ann(exsig) + nbos);
    }

    static constexpr bool contribs_to_frm(size_t exsig, size_t ranksig) {
        return (decode_nfrm_cre(exsig) <= decode_nfrm_cre(ranksig)) && (
                decode_nfrm_cre(ranksig) - decode_nfrm_cre(exsig) == decode_nfrm_ann(ranksig) - decode_nfrm_ann(exsig));
    }

    static constexpr bool contribs_to_bos(size_t exsig, size_t ranksig) {
        return (decode_nbos_cre(exsig) <= decode_nbos_cre(ranksig)) && (
                decode_nbos_cre(ranksig) - decode_nbos_cre(exsig) == decode_nbos_ann(ranksig) - decode_nbos_ann(exsig));
    }

    static constexpr bool contribs_to(size_t exsig, size_t ranksig) {
        return contribs_to_frm(exsig, ranksig) && contribs_to_bos(exsig, ranksig);
    }

    /**
     * @param exsig
     *  excitation signature
     * @return
     *  the exsig representing the hermitian conjugate of the operator represented by the argument
     */
    static constexpr size_t hermconj(size_t exsig) {
        return exsig_utils::encode(decode_nfrm_ann(exsig), decode_nfrm_cre(exsig), decode_nbos_ann(exsig),
                                   decode_nbos_cre(exsig));
    }

    static std::string to_string(size_t exsig) {
        if (exsig > nexsig) return "invalid";
        return std::to_string(decode_nfrm_cre(exsig)) + std::to_string(decode_nfrm_ann(exsig)) +
               std::to_string(decode_nbos_cre(exsig)) + std::to_string(decode_nbos_ann(exsig));
    }

    static constexpr size_t ex_single = encode(1, 1, 0, 0);
    static constexpr size_t ex_double = encode(2, 2, 0, 0);
    static constexpr size_t ex_triple = encode(3, 3, 0, 0);
    static constexpr size_t ex_1100 = ex_single;
    static constexpr size_t ex_2200 = ex_double;
    static constexpr size_t ex_3300 = ex_triple;
    static constexpr size_t ex_1101 = encode(1, 1, 0, 1);
    static constexpr size_t ex_1110 = encode(1, 1, 1, 0);
    static constexpr size_t ex_0001 = encode(0, 0, 0, 1);
    static constexpr size_t ex_0010 = encode(0, 0, 1, 0);
    static constexpr size_t ex_0011 = encode(0, 0, 1, 1);
    static constexpr size_t ex_0022 = encode(0, 0, 2, 2);
}

namespace functor_utils {
    template<typename signature_t, typename T>
    void assert_prototype(){
        static_assert(std::is_convertible<T, std::function<signature_t>>::value,
                      "body function does not conform to required prototype");
    }
    template<typename signature_t, typename T>
    void assert_prototype(const T& t){
        assert_prototype<signature_t, T>();
    }
}


template<typename T>
static std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << utils::to_string(v);
    return os;
}

template<typename T, size_t nind>
static std::ostream &operator<<(std::ostream &os, const std::array<T, nind> &a) {
    os << utils::to_string(array_utils::to_vector(a));
    return os;
}

#endif //M7_UTILS_H
