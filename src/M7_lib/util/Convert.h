//
// Created by rja on 15/06/22.
//

#ifndef M7_CONVERT_H
#define M7_CONVERT_H

#include <map>
#include <iomanip>
#include <stack>
#include <list>
#include "String.h"
#include "Tag.h"

namespace convert {

    constexpr uint_t default_fp(tag::Int<0>) { return 0ul; }

    constexpr uint_t default_fp(tag::Int<1>) { return 6ul; }

    template<typename T>
    constexpr uint_t default_fp() {
        return default_fp(tag::Int<std::is_floating_point<T>::value>());
    }

    /**
     * this should only be delegated to from the below "catch-all" templated overload
     * @tparam T
     *  integer-convertible type which is not explicitly handled by an overload in this namespace. typically an enum
     * @param v
     *  reference to the object being represented as a string
     * @param is_convertible_to_uint_t
     *  static dispatch tag
     * @return
     *  string representation of the integer conversion
     */
    template<typename T>
    static std::string to_string(const T &v, tag::Int<true> /*is_convertible_to_uint_t*/) {
        return std::to_string(static_cast<uint_t>(v));
    }

    /**
     * this should only be delegated to from the below "catch-all" templated overload. typically a user class
     * @tparam T
     *  non integer-convertible type which is not explicitly handled by an overload in this namespace
     * @param v
     *  reference to the object being represented as a string
     * @param is_convertible_to_uint_t
     *  static dispatch tag
     * @return
     *  result of object's own to_string definition
     */
    template<typename T>
    static std::string to_string(const T &v, tag::Int<false> /*is_convertible_to_uint_t*/) {
        return v.to_string();
    }

    template<typename T>
    using string_if_arith_t = typename std::enable_if<std::is_arithmetic<T>::value, std::string>::type;
    template<typename T>
    using string_if_not_arith_t = typename std::enable_if<!std::is_arithmetic<T>::value, std::string>::type;

    /**
     * catch-all for non-arithmetic types: attempt to substitute and call the object's own to_string method, except
     * for types that are convertible to integer
     */
    template<typename T>
    static string_if_not_arith_t<T> to_string(const T &v) {
        return to_string(v, tag::Int<std::is_convertible<T, uint_t>::value>());
    }

    template<typename T>
    static string_if_not_arith_t<T> to_string(const T &v, uint_t /*fp*/) {
        return to_string(v);
    }

    /**
     * catch-all for arithmetic types
     */
    template<typename T>
    static string_if_arith_t<T> to_string(const T &v, uint_t fp = default_fp<T>()) {
        if (v == std::numeric_limits<T>::max()) return "inf";
        std::stringstream tmp;
        tmp << std::scientific << std::setprecision(fp) << v;
        return tmp.str();
    }

    template<typename T>
    static std::string to_string(const T *v) {
        if (!v) return "NULL";
        std::stringstream tmp;
        tmp << v;
        return tmp.str();
    }

    template<typename T>
    static std::string to_string(const T *v, uint_t /*fp*/) {
        return to_string(v);
    }

    template<typename T>
    static std::string to_string(T *const v) {
        if (!v) return "NULL";
        std::stringstream tmp;
        tmp << v;
        return tmp.str();
    }

    template<typename T>
    static std::string to_string(T *const v, uint_t /*fp*/) {
        return to_string(v);
    }

    template<typename T>
    static std::string to_string(const std::complex<T> &v, uint_t fp = default_fp<T>()) {
        return "(" + to_string(v.real(), fp) + ", " + to_string(v.imag(), fp) + ")";
    }

    static std::string to_string(const std::string &str) {
        return "\"" + str + "\"";
    }

    /**
     * such methods must be implemented so we can rely on overloading to work for all types, not just those for which
     * the floating point precision is meaningful
     * @param str
     * @return
     */
    static std::string to_string(const std::string &str, uint_t /*fp*/) {
        return to_string(str);
    }

    template<typename T>
    static std::string to_string(const std::vector<T> &v, uint_t fp = default_fp<T>()) {
        auto fn = [&v, &fp](uint_t i, std::string &word) {
            if (i >= v.size()) return false;
            word = to_string(v[i], fp);
            return true;
        };
        return "[" + string::join(fn, ", ") + "]";
    }

    template<typename key_t, typename value_t>
    static std::string to_string(const std::map<key_t, value_t> &v, uint_t fp = default_fp<value_t>()) {
        auto fn = [&v, &fp](uint_t i, std::string &word) {
            if (i >= v.size()) return false;
            auto it = v.cend();
            std::advance(it, i);
            word = to_string(it->first, fp) + ": " + to_string(it->second, fp);
            return true;
        };
        return "[" + string::join(fn, ", ") + "]";
    }

    template<typename T>
    static std::string to_string(const std::stack<T> &v, uint_t fp = default_fp<T>()) {
        auto cpy = v;
        auto fn = [&cpy, &fp](uint_t i, std::string &word) {
            if (cpy.empty()) return false;
            word = to_string(cpy.top(), fp);
            cpy.pop();
            return true;
        };
        return "[" + string::join(fn, ", ") + "]";
    }

    template<typename T>
    static std::string to_string(const std::list<T> &v, uint_t fp = default_fp<T>()) {
        auto cpy = v;
        std::vector<T> tmp;
        for (const auto &it: v) tmp.push_back(it);
        return to_string(tmp, fp);
    }

    static std::string to_string(bool v) {
        return v ? "true" : "false";
    }

    static std::string to_string(bool v, uint_t) {
        return to_string(v);
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
        for (auto &it: wides) narrows.push_back(convert::safe_narrow<narrow_t>(it));
        return narrows;
    }

    template<typename T1, typename T2>
    static void vector(const std::vector<T1> &v1, std::vector<T2> &v2) {
        static_assert(std::is_convertible<T1, T2>::value, "incompatible types");
        v2.clear();
        v2.reserve(v1.size());
        for (auto &i: v1) v2.push_back(i);
    }
}

template<typename T>
static std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << convert::to_string(v);
    return os;
}


#endif //M7_CONVERT_H