//
// Created by rja on 12/06/22.
//

#include "String.h"

std::string utils::string::join(const std::vector<std::string> &words, const std::string &divider) {
    std::string out;
    for (size_t i = 0ul; i < words.size() - 1; ++i) {
        out += (words[i] + divider);
    }
    out += words[words.size() - 1];
    return out;
}

std::string utils::string::join(const std::vector<std::string> &words) {
    return join(words, " ");
}

std::string utils::string::join(const std::string &word, const size_t &nrepeat, const std::string &divider) {
    return join(std::vector<std::string>(nrepeat, word), divider);
}

std::string utils::string::join(const std::string &word, const size_t &nrepeat) {
    return join(word, nrepeat, " ");
}

std::vector<std::string> utils::string::split(const std::string &line, char delimiter) {
    std::vector<std::string> result{};
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        if (token.size()) result.push_back(token);
    }
    return result;
}

std::vector<std::string> utils::string::split(const std::string &line, const std::string &delimiters) {
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

void utils::string::split(std::string &line, std::vector<std::string> &tokens, const std::string &delimiters) {
    tokens.clear();
    char *ptr;
    ptr = strtok(const_cast<char *>(line.c_str()), delimiters.c_str());
    while (ptr != nullptr) {
        tokens.emplace_back(ptr);
        ptr = strtok(nullptr, delimiters.c_str());
    }
}

std::string utils::string::yn(bool t) {
    return t ? "yes" : "no";
}

std::string utils::string::YN(bool t) {
    return t ? "YES" : "NO";
}

std::string utils::string::memsize(size_t nbyte) {
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

std::string utils::string::boxed(std::string s, size_t padding, char c) {
    std::string res;
    res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
    res += c + std::string(padding, ' ') + s + std::string(padding, ' ') + c + "\n";
    res += std::string(s.size() + 2 * (padding + 1), c) + '\n';
    return res;
}

bool utils::string::is_numeric(const char &c) {
    return '0' <= c && c <= '9';
}

bool utils::string::is_partial_standard_float(const char &c) {
    return is_numeric(c) || c == '.' || c == '-';
}

bool utils::string::is_partial_scientific(const char &c) {
    return is_partial_standard_float(c) || c == 'e' || c == 'E' || c == 'd' || c == 'D' || c == '+';
}

bool utils::string::is_divider(const char &c) {
    return c == ' ' || c == ',' || c == ')' || c == '\r';
}

double utils::string::read_double(const char *&ptr) {
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

size_t utils::string::read_unsigned(const char *&ptr) {
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

int64_t utils::string::read_signed(const char *&ptr) {
    bool pos = true;
    if (*ptr == '-') {
        pos = false;
        ptr++;
    }
    auto tmp = read_unsigned(ptr);
    if (tmp == std::numeric_limits<size_t>::max()) return std::numeric_limits<int64_t>::max();
    return pos ? tmp : -tmp;
}

size_t utils::string::parse_decimal_digit(const char *c) {
    if (*c < '0' || *c > '9') return ~0ul;
    return *c - '0';
}

std::string utils::string::plural(size_t i, std::string plu_ending, std::string sing_ending) {
    return (i == 1) ? sing_ending : plu_ending;
}

std::string utils::string::plural(std::string base, size_t i, std::string plu_ending, std::string sing_ending) {
    return std::to_string(i) + " " + base + plural(i, plu_ending, sing_ending);
}

std::string utils::string::prefix(std::string base, std::string prefix, char delimiter) {
    if (prefix.empty()) return base;
    prefix.push_back(delimiter);
    return prefix+base;
}
