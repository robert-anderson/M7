//
// Created by Robert John Anderson on 2020-01-10.
//

#ifndef M7_UTILS_H
#define M7_UTILS_H

#include "src/defs.h"
#include <vector>
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

    /*
     * tests whether all bytes are zero.
     * This is slow, only use for checking that memsets and other zeroing operations
     * have been successful in debug build
     */
    template<typename T>
    bool is_zero(const T *v) {
        for (size_t ichar = 0ul; ichar < sizeof(T); ++ichar) {
            if (*((char *) v + ichar) != 0) return false;
        }
        return true;
    }

    template<typename T>
    bool is_zero(const T &v) {
        return is_zero(&v);
    }

    template<typename T>
    std::string to_string(const T& v) {
        return std::to_string(v);
    }

    template<typename T>
    std::string to_string(const std::vector<T>& v) {
        std::string string("[");
        for (size_t i = 0ul; i < v.size(); ++i) {
            string += std::to_string(v[i]) + " ";
        }
        string += "]";
        return string;
    }

    template<typename T>
    std::string to_string(const std::stack<T>& v) {
        auto cpy = v;
        std::vector<T> tmp;
        while (!v.empty()){
            tmp.push_back(cpy.top());
            cpy.pop();
        }
        return to_string(tmp);
    }

    template<typename T>
    std::string to_string(const std::list<T>& v) {
        auto cpy = v;
        std::vector<T> tmp;
        for (const auto& i: v) tmp.push_back(i);
        return to_string(tmp);
    }

    template<typename T>
    std::string fp_to_string(const T &v, size_t fp_precision = 6) {
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
    std::string num_to_string(const T &entry, size_t padding = 0, size_t fp_precision = 9) {
        std::string result;
        if (std::is_floating_point<T>::value) result = fp_to_string(entry, fp_precision);
        else if (std::is_integral<T>::value) result = std::to_string(entry);
        return result;
    }

    template<typename T>
    std::string num_to_string(const std::complex<T> &entry, size_t padding = 0, size_t fp_precision = 9) {
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

    template<typename narrow_t, typename wide_t>
    narrow_t safe_narrow(const wide_t &wide) {
        static_assert(std::is_convertible<wide_t, narrow_t>::value, "incompatible types");
        static_assert(sizeof(wide_t) >= sizeof(narrow_t), "wide type must be at least as long as narrow type");
#ifdef SAFE_NARROWING
        if (static_cast<wide_t>(static_cast<narrow_t>(wide)) != wide) throw std::runtime_error("narrowing loses information");
#endif
        return static_cast<narrow_t>(wide);
    }

    template<typename narrow_t, typename wide_t>
    std::vector<narrow_t> safe_narrow(const std::vector<wide_t> &wides) {
        std::vector<narrow_t> narrows;
        narrows.reserve(wides.size());
        for (auto &it : wides) narrows.push_back(utils::safe_narrow<defs::mpi_count>(it));
        return narrows;
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
}

namespace bit_utils {
    template<typename T>
    static inline void clr(T &x, const size_t &i) {
        x &= ~((T) 1ul << i);
    }

    template<typename T>
    static inline void set(T &x, const size_t &i) {
        x |= ((T) 1ul << i);
    }

    template<typename T>
    static inline bool get(const T &x, const size_t &i) {
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

    template<typename T>
    std::string to_string(const T &v){
        std::string tmp;
        for (size_t i=0ul; i<sizeof (T)*CHAR_BIT; ++i) tmp+=get(v, i)?'1':'0';
        return tmp;
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

    static std::vector<std::string> split(const std::string &line, char delimiter) {
        std::vector<std::string> result{};
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, delimiter)) {
            if (token.size()) result.push_back(token);
        }
        return result;
    }

    static std::vector<std::string> split(const std::string &line, const std::string &delimiters) {
        std::string mutable_copy = line;
        std::vector<std::string> result{};
        char *ptr;
        ptr = strtok(const_cast<char *>(mutable_copy.c_str()), delimiters.c_str());
        while (ptr != nullptr) {
            result.emplace_back(ptr);
            ptr = strtok(nullptr, delimiters.c_str());
        }
        return result;
    }

    static std::string yn(bool t) {
        return t ? "yes" : "no";
    }

    static std::string YN(bool t) {
        return t ? "YES" : "NO";
    }

    static std::string memsize(size_t nbyte) {
        if (nbyte < 1e3) {
            return std::to_string(nbyte) + "B";
        } else if (nbyte < 1e6) {
            return std::to_string(nbyte / 1.0e3) + "KB";
        } else if (nbyte < 1e9) {
            return std::to_string(nbyte / (1.0e6)) + "MB";
        } else {
            return std::to_string(nbyte / (1.0e9)) + "GB";
        }
    }

    static std::string boxed(std::string s, size_t padding = 4, char c = '#') {
        std::string res;
        res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
        res += c + std::string(padding, ' ') + s + std::string(padding, ' ') + c + "\n";
        res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
        return res;
    }


    static inline bool is_numeric(const char &c) {
        return '0' <= c && c <= '9';
    }

    static inline bool is_partial_standard_float(const char &c) {
        return is_numeric(c) || c == '.' || c == '-';
    }

    static inline bool is_partial_scientific(const char &c) {
        return is_partial_standard_float(c) || c == 'e' || c == 'E' || c == 'd' || c == 'D' || c == '+';
    }

    static inline bool is_divider(const char &c) {
        return c == ' ' || c == ',' || c == ')' || c=='\r';
    }

    static double read_double(const char *&ptr) {
        const char *begin = nullptr;
        ASSERT(ptr != nullptr)
        for (; *ptr != 0; ptr++) {
            if (!begin) {
                if (is_partial_standard_float(*ptr)) begin = ptr;
            } else {
                if (is_divider(*ptr) && is_numeric(ptr[-1])) {
                    return std::strtod(begin, const_cast<char **>(&ptr));
                } else if (!is_partial_scientific(*ptr)) {
                    begin = nullptr;
                }
            }
        }
        if (begin && is_numeric(ptr[-1])) {
            return std::strtod(begin, const_cast<char **>(&ptr)); // this will decrement the pointer!
        } else {
            return std::numeric_limits<double>::max();
        }
    }

    static size_t read_unsigned(const char *&ptr) {
        const char *begin = nullptr;
        ASSERT(ptr != nullptr)
        for (; *ptr != 0; ptr++) {
            if (!begin) {
                if (is_numeric(*ptr)) begin = ptr;
            } else {
                if (is_divider(*ptr)) {
                    return std::strtoul(begin, const_cast<char **>(&ptr), 10);
                } else if (!is_numeric(*ptr)) {
                    begin = nullptr;
                }
            }
        }
        if (begin && is_numeric(ptr[-1])) {
            return std::strtoul(begin, const_cast<char **>(&ptr), 10); // this will decrement the pointer!
        } else {
            return std::numeric_limits<size_t>::max();
        }
    }

    static int64_t read_signed(const char *&ptr) {
        bool pos = true;
        if (*ptr=='-') {
            pos = false;
            ptr++;
        }
        auto tmp = read_unsigned(ptr);
        if (tmp==std::numeric_limits<size_t>::max()) return std::numeric_limits<int64_t>::max();
        return pos ? tmp : -tmp;
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
}

namespace ci_utils {
    static size_t nalpha(size_t nelec, int spin) {
        size_t spin_odd = std::abs(spin) % 2;
        ASSERT(nelec % 2 == spin_odd)
        size_t nalpha = nelec / 2 + (std::abs(spin)) / 2 + spin_odd;
        return spin >= 0 ? nalpha : nelec - nalpha;
    }

    static size_t nbeta(size_t nelec, int spin) {
        return nelec - nalpha(nelec, spin);
    }

    static size_t fermion_dim(size_t nsite, size_t nelec) {
        return integer_utils::combinatorial(2 * nsite, nelec);
    }

    static size_t fermion_dim(size_t nsite, size_t nelec, int spin) {
        ASSERT(static_cast<size_t>(spin) % 2 == nelec % 2)
        return std::pow(integer_utils::combinatorial(nsite, nalpha(nelec, spin)), 2);
    }

    static size_t boson_dim(size_t boson_nmode, size_t boson_cutoff) {
        return std::pow(boson_cutoff + 1, boson_nmode);
    }
}

namespace mem_utils {

    template<typename T, typename... Args>
    static std::unique_ptr<T> make_unique(Args &&... args) {
        return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
    }

    static void print_cmp(char* c1, char* c2, size_t n){
        for (size_t i=0; i<n; ++i){
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
        }
        else {
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
    for_each(std::tuple<Tp...> &, FuncT&) // Unused arguments are given no names.
    { }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each(std::tuple<Tp...>& t, FuncT& f)
    {
        f(std::get<I>(t));
        for_each<I + 1, FuncT, Tp...>(t, f);
    }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each(const std::tuple<Tp...> &, FuncT&) // Unused arguments are given no names.
    { }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each(const std::tuple<Tp...>& t, FuncT& f)
    {
        f(std::get<I>(t));
        for_each<I + 1, FuncT, Tp...>(t, f);
    }


    /*
     * modifiable / modifiable
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &, std::tuple<Tp...> &, FuncT&) // Unused arguments are given no names.
    { }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...>& t1, std::tuple<Tp...>& t2, FuncT& f)
    {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }


    /*
     * modifiable / const
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...> &, const std::tuple<Tp...> &, FuncT&) // Unused arguments are given no names.
    { }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(std::tuple<Tp...>& t1, const std::tuple<Tp...>& t2, FuncT& f)
    {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }

    /*
     * const / const
     */
    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type
    for_each_pair(const std::tuple<Tp...> &, const std::tuple<Tp...> &, FuncT&) // Unused arguments are given no names.
    { }

    template<std::size_t I = 0, typename FuncT, typename... Tp>
    inline typename std::enable_if<I < sizeof...(Tp), void>::type
    for_each_pair(const std::tuple<Tp...>& t1, const std::tuple<Tp...>& t2, FuncT& f)
    {
        f(std::get<I>(t1), std::get<I>(t2));
        for_each_pair<I + 1, FuncT, Tp...>(t1, t2, f);
    }
}

namespace nd_utils {
    template<typename T>
    T nelement(const std::vector<T>& v){
        T out = 1;
        for (const auto& i: v) out*=i;
        return out;
    };
}

namespace dispatch_utils {
    template<typename T>
    struct TypeTag{};

    template<size_t ind>
    struct IndTag{};

    template<bool t>
    struct BoolTag{};
}

namespace conn_utils {
    static size_t left(const size_t& ispinorb, const size_t& nsite) {
        if (ispinorb==0 || ispinorb==nsite) return ~0ul;
        return ispinorb-1;
    }
    static size_t left_pbc(const size_t& ispinorb, const size_t& nsite) {
        if (ispinorb==0) return nsite-1;
        else if (ispinorb==nsite) return 2*nsite-1;
        return ispinorb-1;
    }
    static size_t right(const size_t& ispinorb, const size_t& nsite) {
        if (ispinorb==nsite-1 || ispinorb==2*nsite-1) return ~0ul;
        return ispinorb+1;
    }
    static size_t right_pbc(const size_t& ispinorb, const size_t& nsite) {
        if (ispinorb==nsite-1) return 0;
        else if (ispinorb==2*nsite-1) return nsite;
        return ispinorb+1;
    }
}

namespace foreach {
    /**
     * "compile time number of dimensions"
     */
    namespace ctnd {
        template<size_t nind>
        using inds_t = std::array<size_t, nind>;
        template<size_t nind>
        using body_fn_t = std::function<void(const inds_t<nind> &)>;
        template<size_t nind>
        using cr_body_fn_t = const body_fn_t<nind> &;

        template<size_t nind>
        struct Base {
            const size_t m_nterm;
            inds_t<nind> m_inds;

            Base(size_t nterm) : m_nterm(nterm) {}
        };

        template<size_t nind>
        struct Unrestricted : Base<nind> {
            const inds_t<nind> m_shape;
        private:
            static size_t nterm(const inds_t<nind> &shape) {
                if (!nind) return 0;
                size_t n = 1;
                for (size_t i = 0ul; i < nind; ++i) n *= shape[i];
                return n;
            }

        public:
            Unrestricted(const inds_t<nind> &shape) : Base<nind>(nterm(shape)), m_shape(shape) {}

            template<size_t ilevel>
            void level_loop(cr_body_fn_t<nind> body, dispatch_utils::IndTag<ilevel>) {
                constexpr size_t iind = ilevel - 1;
                auto &ind = Base<nind>::m_inds[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) level_loop(body, dispatch_utils::IndTag<ilevel + 1>());
            }

            void level_loop(cr_body_fn_t<nind> body, dispatch_utils::IndTag<nind>) {
                constexpr size_t iind = nind - 1;
                auto &ind = Base<nind>::m_inds[iind];
                const auto &extent = m_shape[iind];
                for (ind = 0ul; ind < extent; ++ind) body(Base<nind>::m_inds);
            }

            void loop(cr_body_fn_t<nind> body, dispatch_utils::BoolTag<true>) {}

            void loop(cr_body_fn_t<nind> body, dispatch_utils::BoolTag<false>) {
                level_loop(body, dispatch_utils::IndTag<1>());
            }

            void operator()(cr_body_fn_t<nind> body) {
                loop(body, dispatch_utils::BoolTag<nind == 0>());
            }
        };

        template<size_t nind, bool strict = true, bool ascending = true>
        struct Ordered : Base<nind> {
            using Base<nind>::m_inds;
            using inds_t = std::array<size_t, nind>;
            const size_t m_n;

        public:
            Ordered(const size_t &n) : Base<nind>(100), m_n(n) {}

            template<size_t ilevel>
            void level_loop(cr_body_fn_t<nind> body, dispatch_utils::IndTag<ilevel>) {
                constexpr size_t iind = ascending ? (nind - ilevel) : (ilevel - 1);
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = (iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1]) + !strict;
                for (ind = 0ul; ind < extent; ++ind) level_loop(body, dispatch_utils::IndTag<ilevel + 1>());
            }

            void level_loop(cr_body_fn_t<nind> body, dispatch_utils::IndTag<nind>) {
                constexpr size_t iind = ascending ? 0 : nind - 1;
                constexpr size_t iind_unrestrict = ascending ? nind - 1 : 0;
                auto &ind = Base<nind>::m_inds[iind];
                const auto extent = (iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1]) + !strict;
                for (ind = 0ul; ind < extent; ++ind) body(m_inds);
            }

            void loop(cr_body_fn_t<nind> body, dispatch_utils::BoolTag<true>) {}

            void loop(cr_body_fn_t<nind> body, dispatch_utils::BoolTag<false>) {
                level_loop(body, dispatch_utils::IndTag<1>());
            }

            void operator()(cr_body_fn_t<nind> body) {
                loop(body, dispatch_utils::BoolTag<nind == 0>());
            }
        };
    }

    /**
     * "run time number of dimensions"
     */
    namespace rtnd {
        using inds_t = std::vector<size_t>;
        using body_fn_t = std::function<void(const inds_t &)>;
        using cr_body_fn_t = const body_fn_t&;

        struct Base {
            const size_t m_nind;
            const size_t m_nterm;
            inds_t m_inds;

            Base(size_t nind, size_t nterm) : m_nind(nind),
                                              m_nterm(nterm), m_inds(m_nind, 0) {}
        };

        struct Unrestricted : Base {
            const inds_t m_shape;
        private:
            static size_t nterm(const inds_t &shape) {
                if (shape.empty()) return 0;
                size_t n = 1;
                for (const auto &i: shape) n *= i;
                return n;
            }

        public:
            Unrestricted(const inds_t &shape) : Base(shape.size(), nterm(shape)), m_shape(shape) {}

            void level_loop(cr_body_fn_t body, size_t ilevel) {
                const auto &iind = ilevel - 1;
                auto &ind = m_inds[iind];
                const auto &extent = m_shape[iind];
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(body, ilevel + 1);
                else
                    for (ind = 0ul; ind < extent; ++ind) body(m_inds);
            }

            void operator()(cr_body_fn_t body) {
                if (m_nind == 0) return;
                level_loop(body, 1);
            }
        };

        template<bool strict = true, bool ascending = true>
        struct Ordered : Base {
            const size_t m_n;
        public:
            Ordered(const size_t &n, const size_t &r) :
                Base(r, integer_utils::combinatorial(strict ? n : (n+r)-1, r)), m_n(n) {}

            void level_loop(cr_body_fn_t body, size_t ilevel) {
                const size_t iind = ascending ? (m_nind - ilevel) : (ilevel - 1);
                const size_t iind_unrestrict = ascending ? m_nind - 1 : 0;
                auto &ind = m_inds[iind];
                const auto extent = (iind == iind_unrestrict ? m_n : m_inds[ascending ? iind + 1 : iind - 1]) + !strict;
                if (ilevel < m_nind)
                    for (ind = 0ul; ind < extent; ++ind) level_loop(body, ilevel + 1);
                else
                    for (ind = 0ul; ind < extent; ++ind) body(m_inds);
            }

            void operator()(cr_body_fn_t body) {
                if (m_nind == 0) return;
                level_loop(body, 1);
            }
        };


    }
}

#endif //M7_UTILS_H
