//
// Created by rja on 12/06/22.
//

#include <cstring>
#include "String.h"

str_t string::join(const strv_t &words, const str_t &divider) {
    auto fn = [&words](uint_t i, str_t& word) {
        if (i >= words.size()) return false;
        word = words[i];
        return true;
    };
    return join(fn, divider);
}

str_t string::join(const strv_t &words) {
    return join(words, " ");
}

str_t string::join(const str_t &word, const uint_t &nrepeat, const str_t &divider) {
    return join(strv_t(nrepeat, word), divider);
}

str_t string::join(const str_t &word, const uint_t &nrepeat) {
    return join(word, nrepeat, " ");
}

strv_t string::split(const str_t &line, char delimiter) {
    strv_t result{};
    std::stringstream ss(line);
    str_t token;
    while (std::getline(ss, token, delimiter)) {
        if (token.size()) result.push_back(token);
    }
    return result;
}

strv_t string::split(const str_t &line, const str_t &delimiters) {
    str_t mutable_copy = line;
    strv_t result{};
    char *ptr;
    ptr = std::strtok(const_cast<char *>(mutable_copy.c_str()), delimiters.c_str());
    while (ptr != nullptr) {
        result.emplace_back(ptr);
        ptr = std::strtok(nullptr, delimiters.c_str());
    }
    return result;
}

void string::split(str_t &line, strv_t &tokens, const str_t &delimiters) {
    tokens.clear();
    char *ptr;
    ptr = std::strtok(const_cast<char *>(line.c_str()), delimiters.c_str());
    while (ptr != nullptr) {
        tokens.emplace_back(ptr);
        ptr = std::strtok(nullptr, delimiters.c_str());
    }
}

str_t string::yn(bool t) {
    return t ? "yes" : "no";
}

str_t string::YN(bool t) {
    return t ? "YES" : "NO";
}

str_t string::memsize(uint_t nbyte) {
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

str_t string::boxed(str_t s, uint_t padding, char c) {
    str_t res;
    res += str_t(s.size() + 2 * (padding + 1), c) + '\n';
    res += c + str_t(padding, ' ') + s + str_t(padding, ' ') + c + "\n";
    res += str_t(s.size() + 2 * (padding + 1), c) + '\n';
    return res;
}

bool string::is_numeric(const char &c) {
    return '0' <= c && c <= '9';
}

bool string::is_partial_standard_float(const char &c) {
    return is_numeric(c) || c == '.' || c == '-';
}

bool string::is_partial_scientific(const char &c) {
    return is_partial_standard_float(c) || c == 'e' || c == 'E' || c == 'd' || c == 'D' || c == '+';
}

bool string::is_divider(const char &c) {
    return c == ' ' || c == ',' || c == ')' || c == '\r';
}

double string::read_double(const char *&ptr) {
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

uint_t string::read_unsigned(const char *&ptr) {
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
        return std::numeric_limits<uint_t>::max();
    }
}

int64_t string::read_signed(const char *&ptr) {
    bool pos = true;
    if (*ptr == '-') {
        pos = false;
        ptr++;
    }
    auto tmp = read_unsigned(ptr);
    if (tmp == std::numeric_limits<uint_t>::max()) return std::numeric_limits<int64_t>::max();
    return pos ? tmp : -tmp;
}

uint_t string::parse_decimal_digit(const char *c) {
    if (*c < '0' || *c > '9') return ~0ul;
    return *c - '0';
}

str_t string::plural(uint_t i, str_t plu_ending, str_t sing_ending) {
    return (i == 1) ? sing_ending : plu_ending;
}

str_t string::plural(str_t base, uint_t i, str_t plu_ending, str_t sing_ending) {
    return std::to_string(i) + " " + base + plural(i, plu_ending, sing_ending);
}

str_t string::prefix(str_t base, str_t prefix, char delimiter) {
    if (prefix.empty()) return base;
    prefix.push_back(delimiter);
    return prefix+base;
}

void string::to_upper(str_t &str) {
    const auto shift = 'A' - 'a';
    for (auto& c: str) if (c >= 'a' && c <= 'z') c += shift;
}

void string::to_lower(str_t &str) {
    const auto shift = 'a' - 'A';
    for (auto& c: str) if (c >= 'A' && c <= 'Z') c += shift;
}
