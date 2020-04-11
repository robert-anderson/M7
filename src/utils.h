//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_UTILS_H
#define M7_UTILS_H

#include <vector>
#include <complex>
#include <iostream>
#include <cmath>
#include <assert.h>
#include <limits>
#include <numeric>
#include <x86intrin.h>
#include <climits>
#include <iomanip>

namespace utils {
    template <typename T>
    std::string to_string(T v, size_t n){
        std::string string("[");
        for (size_t i=0ul; i<n; ++i){
            string+=std::to_string(v[i])+" ";
        }
        string+="]";
        return string;
    }

    template <typename T>
    std::string to_string(T v) {
        return to_string(v, v.size());
    }

    template<typename T>
    std::string fp_to_string(const T &v, size_t fp_precision=6){
        assert(std::is_floating_point<T>::value);
        std::stringstream tmp;
        tmp << std::scientific << std::setprecision(fp_precision) << v;
        return tmp.str();
    }

    template<typename T>
    std::string num_to_string(const T &entry, size_t padding = 0, size_t fp_precision=6) {
        std::string result;
        if (std::is_floating_point<T>::value) result = fp_to_string(entry, fp_precision);
        else if (std::is_integral<T>::value) result = std::to_string(entry);
        auto decimal_length = std::numeric_limits<T>::digits10;
        //assert(result.size()<=decimal_length);
        //result.insert(result.begin(), padding + decimal_length - result.size(), ' ');
        return result;
    }

    template<typename T>
    std::string num_to_string(const std::complex<T> &entry, size_t padding = 0, size_t fp_precision=6) {
        auto tmp_string = fp_to_string(entry.real(), fp_precision) +
                          (entry.imag() < 0 ? "" : "+") + fp_to_string(entry.imag(), fp_precision) + "i";
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }


    template <typename T>
    void print(T v){
        std::cout << to_string(v) << std::endl;
    }

    template <typename T>
    void print(T v, size_t n){
        std::cout << to_string(v, n) << std::endl;
    }
}

namespace integer_utils {

    template<typename T>
    static typename std::enable_if<std::is_integral<T>::value, T>::type
    divceil(const T& num, const T& denom){
        return num%denom ? num/denom+1:num/denom;
    }

    static size_t factorial(const size_t &n){
        assert(n<((size_t)-1)/2);
        size_t out = 1ul;
        if (n<1) return 1ul;
        for (size_t i = 1ul; i<=n; ++i) out*=i;
        return out;
    }

    static size_t combinatorial(const size_t &n, const size_t &r){
        /*
         * n choose r = n! / ((n-r)!r!)
         * compute numerator and denominator simultaneously whenever an
         * exact quotient can be computed to avoid premature overflow
         */
        assert(n>=r);
        if (r == 0) return 1ul;
        if (n == 1) return 1ul;
        if (r == n) return 1ul;

        size_t out = 1ul;
        size_t ni = 0ul;
        size_t ri = 0ul;
        while (1){
            if (ri<r && out%(r-ri)==0) {
                out/=r-(ri++);
            }
            else out*=n-(ni++);
            assert(ni<=r); // overflow occurred.
            if (ri==r && ni==r) return out;
        }
    }
}

namespace bit_utils{
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

    template <typename T>
    static inline size_t next_setbit(T& work);

    template <>
    inline size_t next_setbit(unsigned long long &work) {
        static_assert(sizeof(typeof(work))==8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    template <>
    inline size_t next_setbit(unsigned long &work) {
        static_assert(sizeof(typeof(work))==8, "Data length not supported");
        size_t result = __tzcnt_u64(work);
        bit_utils::clr(work, result);
        return result;
    }

    template <>
    inline size_t next_setbit(unsigned &work) {
        static_assert(sizeof(typeof(work))==4, "Data length not supported");
        size_t result = __tzcnt_u32(work);
        bit_utils::clr(work, result);
        return result;
    }

    template <typename T>
    static inline size_t nsetbit(const T& work);

    template <>
    inline size_t nsetbit(const unsigned long long& work){
        static_assert(sizeof(typeof(work))==8, "Data length not supported");
        return _popcnt64(work);
    }

    template <>
    inline size_t nsetbit(const unsigned long& work){
        static_assert(sizeof(typeof(work))==8, "Data length not supported");
        return _popcnt64(work);
    }

    template <>
    inline size_t nsetbit(const unsigned& work){
        static_assert(sizeof(typeof(work))==4, "Data length not supported");
        return _popcnt32(work);
    }

    template <typename T>
    T truncate(T &v, size_t n){
        const auto nbit = sizeof(T)*CHAR_BIT-n;
        return (v<<nbit)>>nbit;
    }

}


namespace string_utils {
    static std::string join(const std::vector<std::string> &words, const std::string &divider, const bool &bookends) {
        std::string out{""};
        if (bookends) out += divider;
        for (size_t i=0ul; i < words.size() - 1; ++i) {
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
    void normalize(std::vector<T> &v, T norm=T(1)){
        T tot = std::accumulate(v.begin(), v.end(), T(0));
        T fac = norm/tot;
        for (auto &i:v) i*=fac;
    }
}

#endif //M7_UTILS_H
