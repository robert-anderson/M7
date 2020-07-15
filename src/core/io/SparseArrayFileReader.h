//
// Created by rja on 15/07/2020.
//

#ifndef M7_SPARSEARRAYFILEREADER_H
#define M7_SPARSEARRAYFILEREADER_H

#include <src/core/util/consts.h>
#include <iostream>
#include <src/core/util/defs.h>
#include <src/core/util/utils.h>
#include <algorithm>
#include "FileReader.h"

template<typename T>
class SparseArrayFileReader : public FileReader {
    const size_t m_nind;
    bool m_indsfirst;
    bool m_complex_valued;
public:
    SparseArrayFileReader(const std::string &fname, size_t nind) :
            FileReader(fname), m_nind(nind) {
        size_t iline = first_valid_line(fname, m_nind, m_indsfirst, m_complex_valued);
        if (iline == ~0ul) throw std::runtime_error("No valid entries found");
        if (m_complex_valued && !consts::is_complex<T>())
            throw std::runtime_error("Trying to read complex-valued array entries into a real container");
        else if (consts::is_complex<T>() && !m_complex_valued)
            std::cout << "Reading real-valued array into complex container, consider recompiling with real arithmetic"
                      << std::endl;
        skip(iline);
    }

    bool next(defs::inds &inds, T &v) {
        std::string line;
        bool result = FileReader::next(line);
        bool complex_valued = line.find(',') != ~0ul;
        if (!result) return false;
        extract(line, m_nind, m_indsfirst, complex_valued, inds, v);
        return true;
    }

    size_t first_valid_line(const std::string &fname, const size_t &nind, bool &indsfirst, bool &complex_valued) {
        FileReader iterator(fname);
        std::string line;
        while (iterator.next(line)) {
            if (validate(line, nind, indsfirst, complex_valued)) return iterator.iline();
        }
        return ~0ul; // no valid line found
    }

    static bool validate(const std::string &line, const size_t &nind, bool &indsfirst, bool &complex_valued) {
        // can only be complex valued if the real and imag parts are delimited by a comma
        complex_valued = line.find(',')!=~0ul;
        auto tokens = string_utils::split(line, "() ,");
        if (tokens.size() > nind + 1 + complex_valued) return false;
        std::vector<double> doubles{};
        for (auto token : tokens) {
            try {
                doubles.push_back(std::stod(token));
            }
            catch (std::invalid_argument &) {
                return false;
            }
        }
        using namespace float_utils;
        if (!is_integral(doubles[0]) || (complex_valued && !is_integral(doubles[1]))) {
            // float(s) come first
            indsfirst = false;
            auto begin = doubles.data() + (complex_valued ? 2 : 1);
            auto end = begin + nind;
            if (!std::all_of(begin, end, is_integral<double>)) return false;
        } else if (!is_integral(doubles[doubles.size() - 1]) ||
                   (complex_valued && !is_integral(doubles[doubles.size() - 2]))) {
            // float(s) come last
            indsfirst = true;
            auto begin = doubles.data();
            auto end = begin + nind;
            if (!std::all_of(begin, end, is_integral<double>)) return false;
        }
        return true;
    }

    static void
    extract(const std::string &line, size_t nind, bool indsfirst, bool complex_valued, defs::inds &inds, T &value) {
        const char *p = line.begin().base();
        auto f = [&p, &nind, &complex_valued, &inds, &value](bool indsfirst) {
            if (indsfirst) {
                for (size_t i = 0ul; i < nind; ++i) {
                    inds[i] = read_unsigned(p);
                }
            } else {
                value = read_double(p);
                if (consts::is_complex<T>() && complex_valued)
                    complex_utils::set_imag_part(value, read_double(p));
            }
        };
        f(indsfirst);
        f(!indsfirst);
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
        return c == ' ' || c == ',' || c == ')';
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
};



#endif //M7_SPARSEARRAYFILEREADER_H
