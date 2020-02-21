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

namespace utils {
    template <typename T>
    std::string to_string(T v){
        std::string string("[");
        for (auto i : v){
            string+=std::to_string(i)+" ";
        }
        string+="]";
        return string;
    }
    template<typename T>
    std::string num_to_string(T &entry, size_t padding = 0) {
        auto tmp_string = std::to_string(entry);
        auto decimal_length = std::numeric_limits<T>::digits10;
        tmp_string.insert(tmp_string.begin(), padding + decimal_length - tmp_string.size(), ' ');
        return tmp_string;
    }

    template<typename T>
    std::string num_to_string(std::complex<T> &entry, size_t padding = 0) {
        auto tmp_string = std::to_string(entry.real()) +
                          (entry.imag() < 0 ? "" : "+") + std::to_string(entry.imag()) + "i";
        tmp_string.insert(tmp_string.begin(), padding, ' ');
        return tmp_string;
    }


    template <typename T>
    void print(T v){
        std::cout << to_string(v) << std::endl;
    }
}

namespace integer_utils {
    static size_t factorial(const size_t &n){
        assert(n<((size_t)-1)/2);
        size_t out{1ul};
        if (n<1) return 1ul;
        for (auto i{1ul}; i<=n; ++i) out*=i;
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

        size_t out{1ul};
        auto ni{0ul};
        auto ri{0ul};
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


namespace string_utils {
    static std::string join(const std::vector<std::string> &words, const std::string &divider, const bool &bookends) {
        std::string out{""};
        if (bookends) out += divider;
        for (auto i{0ul}; i < words.size() - 1; ++i) {
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


#endif //M7_UTILS_H
