//
// Created by rja on 12/06/22.
//

#ifndef M7_UTIL_STRING_H
#define M7_UTIL_STRING_H

#include "utils.h"

namespace utils {
    namespace string {
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
}


#endif //M7_UTIL_STRING_H
