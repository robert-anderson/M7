//
// Created by rja on 19/05/2020.
//

#include "utils.h"


size_t integer_utils::rectmap(const size_t &irow, const size_t &icol, const size_t &ncol) {
    // rectangular map
    /*
     * n         i,j
     * -----------------
     * 0 1 2 3   0,0 0,1 0,2 0,3
     * 4 5 6 7   1,0 1,1 1,2 1,3
     */
    return irow * ncol + icol;
}

void integer_utils::inv_rectmap(size_t &irow, size_t &icol, const size_t &ncol, const size_t &flat) {
    irow = flat / ncol;
    icol = flat - irow * ncol;
}

size_t integer_utils::trigmap(const size_t &i, const size_t &j) {
    ASSERT(i >= j);
    /*
     * n        i,j
     * -----------------
     * 0        0,0
     * 1 2      1,0 1,1
     * 3 4 5    2,0 2,1 2,2
     * 6 7 8 9  3,0 3,1 3,2 3,3
     */
    return (i * (i + 1)) / 2 + j;
}

size_t integer_utils::npair(const size_t &ndim) {
    return trigmap(ndim, 0);
}

void integer_utils::inv_trigmap(size_t &i, size_t &j, const size_t &n) {
    i = (size_t) ((std::sqrt(1 + 8 * (double) n) - 1) / 2);
    j = n - (i * (i + 1)) / 2;
}

size_t integer_utils::strigmap(const size_t &i, const size_t &j) {
    ASSERT(i > j);
    // strict triangular map i>j
    /*
     * n        i,j
     * -----------------
     * 0        1,0
     * 1 2      2,0 2,1
     * 3 4 5    3,0 3,1 3,2
     * 6 7 8 9  4,0 4,1 4,2 4,3
     */
    return (i * (i - 1)) / 2 + j;
}

void integer_utils::inv_strigmap(size_t &i, size_t &j, const size_t &n) {
    i = (size_t) ((std::sqrt(1 + 8 * (double) n) + 1) / 2);
    j = n - (i * (i - 1)) / 2;
}

size_t integer_utils::nspair(const size_t &ndim) {
    if (!ndim) return 0ul;
    return strigmap(ndim, 0);
}

size_t integer_utils::factorial(const size_t &n) {
    ASSERT(n < ((size_t) -1) / 2);
    size_t out = 1ul;
    if (n < 1) return 1ul;
    for (size_t i = 1ul; i <= n; ++i) out *= i;
    return out;
}

size_t integer_utils::combinatorial(const size_t &n, const size_t &r) {
    /*
     * n choose r = n! / ((n-r)!r!)
     * compute numerator and denominator simultaneously whenever an
     * exact quotient can be computed to avoid premature overflow
     */
    ASSERT(n >= r);
    if (r == 0) return 1ul;
    if (n == 1) return 1ul;
    if (r == n) return 1ul;

    size_t out = 1ul;
    size_t ni = 0ul;
    size_t ri = 0ul;
    while (1) {
        if (ri < r && out % (r - ri) == 0) {
            out /= r - (ri++);
        } else out *= n - (ni++);
        ASSERT(ni <= r); // overflow occurred.
        if (ri == r && ni == r) return out;
    }
}

size_t integer_utils::combinatorial_with_repetition(const size_t &n, const size_t &r) {
    return combinatorial(n + r - 1, r);
}

std::string string_utils::join(const std::vector<std::string> &words, const std::string &divider) {
    std::string out;
    for (size_t i = 0ul; i < words.size() - 1; ++i) {
        out += (words[i] + divider);
    }
    out += words[words.size() - 1];
    return out;
}

std::string string_utils::join(const std::vector<std::string> &words) {
    return join(words, " ");
}

std::string string_utils::join(const std::string &word, const size_t &nrepeat, const std::string &divider) {
    return join(std::vector<std::string>(nrepeat, word), divider);
}

std::string string_utils::join(const std::string &word, const size_t &nrepeat) {
    return join(word, nrepeat, " ");
}

std::vector<std::string> string_utils::split(const std::string &line, char delimiter) {
    std::vector<std::string> result{};
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        if (token.size()) result.push_back(token);
    }
    return result;
}

std::vector<std::string> string_utils::split(const std::string &line, const std::string &delimiters) {
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

void string_utils::split(std::string &line, std::vector<std::string> &tokens, const std::string &delimiters) {
    tokens.clear();
    char *ptr;
    ptr = strtok(const_cast<char *>(line.c_str()), delimiters.c_str());
    while (ptr != nullptr) {
        tokens.emplace_back(ptr);
        ptr = strtok(nullptr, delimiters.c_str());
    }
}

std::string string_utils::yn(bool t) {
    return t ? "yes" : "no";
}

std::string string_utils::YN(bool t) {
    return t ? "YES" : "NO";
}

std::string string_utils::memsize(size_t nbyte) {
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

std::string string_utils::boxed(std::string s, size_t padding, char c) {
    std::string res;
    res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
    res += c + std::string(padding, ' ') + s + std::string(padding, ' ') + c + "\n";
    res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
    return res;
}

bool string_utils::is_numeric(const char &c) {
    return '0' <= c && c <= '9';
}

bool string_utils::is_partial_standard_float(const char &c) {
    return is_numeric(c) || c == '.' || c == '-';
}

bool string_utils::is_partial_scientific(const char &c) {
    return is_partial_standard_float(c) || c == 'e' || c == 'E' || c == 'd' || c == 'D' || c == '+';
}

bool string_utils::is_divider(const char &c) {
    return c == ' ' || c == ',' || c == ')' || c == '\r';
}

double string_utils::read_double(const char *&ptr) {
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

size_t string_utils::read_unsigned(const char *&ptr) {
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

int64_t string_utils::read_signed(const char *&ptr) {
    bool pos = true;
    if (*ptr == '-') {
        pos = false;
        ptr++;
    }
    auto tmp = read_unsigned(ptr);
    if (tmp == std::numeric_limits<size_t>::max()) return std::numeric_limits<int64_t>::max();
    return pos ? tmp : -tmp;
}

size_t string_utils::parse_decimal_digit(const char *c) {
    if (*c < '0' || *c > '9') return ~0ul;
    return *c - '0';
}

std::string string_utils::plural(size_t i, std::string plu_ending, std::string sing_ending) {
    return (i == 1) ? sing_ending : plu_ending;
}

std::string string_utils::plural(std::string base, size_t i, std::string plu_ending, std::string sing_ending) {
    return std::to_string(i) + " " + base + plural(i, plu_ending, sing_ending);
}
