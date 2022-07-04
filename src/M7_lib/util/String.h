//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_STRING_H
#define M7_UTIL_STRING_H

#include "Functor.h"

namespace string {

    template<typename fn_t>
    str_t join(const fn_t &word_getter_fn, const str_t &divider) {
        functor::assert_prototype<bool(uint_t, str_t &)>(word_getter_fn);
        str_t out;
        str_t tmp;
        if (!word_getter_fn(0ul, tmp)) return {};
        out.insert(out.end(), tmp.cbegin(), tmp.cend());
        for (uint_t i = 1ul;; ++i) {
            if (!word_getter_fn(i, tmp)) return out;
            out.insert(out.end(), divider.cbegin(), divider.cend());
            out.insert(out.end(), tmp.cbegin(), tmp.cend());
        }
        // shouldn't reach here
    }

    str_t join(const strv_t &words, const str_t &divider);

    str_t join(const strv_t &words);

    str_t join(const str_t &word, const uint_t &nrepeat, const str_t &divider);

    str_t join(const str_t &word, const uint_t &nrepeat);

    strv_t split(const str_t &line, char delimiter);

    strv_t split(const str_t &line, const str_t &delimiters);

    void split(str_t &line, strv_t &tokens, const str_t &delimiters);

    str_t yn(bool t);

    str_t YN(bool t);

    str_t memsize(uint_t nbyte);

    str_t boxed(str_t s, uint_t padding = 4, char c = '#');


    inline bool is_numeric(const char &c);

    inline bool is_partial_standard_float(const char &c);

    inline bool is_partial_scientific(const char &c);

    inline bool is_divider(const char &c);

    double read_double(const char *&ptr);

    uint_t read_unsigned(const char *&ptr);

    int64_t read_signed(const char *&ptr);

    uint_t parse_decimal_digit(const char *c);

    str_t plural(uint_t i, str_t plu_ending = "s", str_t sing_ending = "");

    str_t plural(str_t base, uint_t i, str_t plu_ending = "s", str_t sing_ending = "");

    str_t prefix(str_t base, str_t prefix = "", char delimiter = ' ');

}

#endif //M7_UTIL_STRING_H
