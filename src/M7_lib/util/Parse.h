//
// Created by anderson on 20/02/2023.
//

#ifndef M7_PARSE_H
#define M7_PARSE_H

#include "M7_lib/io/Logging.h"
#include "Arith.h"

/**
 * parsing is done in four modes:
 *  throwing: throws std::invalid_argument and std::out_of_range exceptions if parse fails
 *  catching: catches the above exceptions and returns false in the case that an instance of either was raised
 *  debug_checked: raises an error in debug mode if instance of either exception is thrown
 *  checked: raises an error if instance of either exception is thrown
 */
namespace parse {
    static void throwing(const str_t& str, uint_t& v) {
        v = std::stoul(str);
    }

    static void throwing(const str_t& str, long& v) {
        v = std::stol(str);
    }

    static void throwing(const str_t& str, int& v) {
        v = std::stoi(str);
    }

    static void throwing(const str_t& str, double& v) {
        v = std::stod(str);
    }

    typedef strv_t::const_iterator c_iter_token_t;

    /**
     * non-complex number from one string
     */
    template<typename T>
    static void throwing(c_iter_token_t begin, c_iter_token_t /*end*/, T& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        return throwing(*begin, v);
    }

    /**
     * complex number from either one string or pair
     */
    template<typename T>
    static void throwing(c_iter_token_t begin, c_iter_token_t end, std::complex<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        with_except(*begin, arith::real_ref(v));
        if (begin+1==end) return;
        with_except(*(begin+1), arith::imag_ref(v));
    }

    template<typename T>
    static void throwing(c_iter_token_t begin, c_iter_token_t end, v_t<T>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it!=end; ++it) {
            v.emplace_back();
            throwing(it, it+1, v.back());
        }
    }

    template<typename T>
    static void throwing(c_iter_token_t begin, c_iter_token_t end, v_t<std::complex<T>>& v){
        static_assert(std::is_arithmetic<T>::value, "can only parse arithmetic types");
        v.clear();
        for (auto it = begin; it != end; ++it) {
            v.emplace_back();
            throwing(*it, arith::real_ref(v.back()));
            if (++it == end) return;
            throwing(*it, arith::imag_ref(v.back()));
        }
    }

    /*
     * parse many from many strings
     */
    template<typename T>
    void throwing(const strv_t& strs, v_t<T>& v) {
        if (v.size() < strs.size()) v.resize(strs.size());
        uint_t i = 0ul;
        for (auto& str: strs) throwing(str, v[i++]);
    }

    template<typename T>
    bool catching(const str_t& str, T& v) {
        try {
            throwing(str, v);
        }
        catch (const std::invalid_argument& ex){
            return false;
        }
        catch (const std::out_of_range& ex){
            return false;
        }
        return true;
    }


    template<typename T>
    bool catching(c_iter_token_t begin, c_iter_token_t end, T& v) {
        try {
            throwing(begin, end, v);
        }
        catch (const std::invalid_argument& ex){
            return false;
        }
        catch (const std::out_of_range& ex){
            return false;
        }
        return true;
    }

    template<typename T>
    bool catching(const strv_t& strs, v_t<T>& v) {
        try {
            throwing(strs, v);
        }
        catch (const std::invalid_argument& ex){
            return false;
        }
        catch (const std::out_of_range& ex){
            return false;
        }
        return true;
    }

    template<typename T>
    void debug_checked(const str_t& str, T& v) {
        const auto success = catching(str, v);
        DEBUG_ONLY(success);
        DEBUG_ASSERT_TRUE(success, logging::format("string \"{}\" could not be parsed", str));
    }

    template<typename T>
    void debug_checked(c_iter_token_t begin, c_iter_token_t end, T& v) {
        const auto success = catching(begin, end, v);
        DEBUG_ONLY(success);
        DEBUG_ASSERT_TRUE(success, logging::format(
                "strings \"{}\", \"{}\" could not be parsed", *begin, *end));
    }

    template<typename T>
    void debug_checked(const strv_t& strs, v_t<T>& v) {
        const auto success = catching(strs, v);
        DEBUG_ONLY(success);
        DEBUG_ASSERT_TRUE(success, logging::format("strings \"{}\" could not be parsed", strs));
    }

    template<typename T>
    void checked(const str_t& str, T& v) {
        const auto success = catching(str, v);
        REQUIRE_TRUE(success, logging::format("string \"{}\" could not be parsed", str));
    }

    template<typename T>
    void checked(c_iter_token_t begin, c_iter_token_t end, T& v) {
        const auto success = catching(begin, end, v);
        REQUIRE_TRUE(success, logging::format(
                "strings \"{}\", \"{}\" could not be parsed", *begin, *end));
    }

    template<typename T>
    void checked(const strv_t& strs, v_t<T>& v) {
        const auto success = catching(strs, v);
        REQUIRE_TRUE(success, logging::format("strings \"{}\" could not be parsed", strs));
    }
}


#endif //M7_PARSE_H
