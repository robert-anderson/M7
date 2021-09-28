//
// Created by rja on 15/07/2020.
//

#ifndef M7_SPARSEARRAYFILEREADER_H
#define M7_SPARSEARRAYFILEREADER_H

#include <src/core/util/consts.h>
#include <iostream>
#include <src/defs.h>
#include <src/core/util/utils.h>
#include <algorithm>
#include <src/core/util/Ternary.h>
#include "FileReader.h"
#include "Logging.h"
#include "src/core/parallel/MPIAssert.h"


template<typename T>
class SparseArrayFileReader : public FileReader {
    const size_t m_nind;
public:
    Tern m_indsfirst;
    Tern m_complex_valued;

    SparseArrayFileReader(const std::string &fname, size_t nind,
                          Tern indsfirst = Tern(), Tern complex_valued = Tern()) :
            FileReader(fname), m_nind(nind), m_indsfirst(indsfirst), m_complex_valued(complex_valued) {
        reset();
        if (m_indsfirst) log::info("Reading sparse array from file with indices before values");
        else log::info("Reading sparse array from file with indices after values");
        REQUIRE_TRUE_ALL(!m_complex_valued || consts::is_complex<T>(), "Trying to read complex-valued array entries into a real container");
        if (consts::is_complex<T>() && !m_complex_valued)
            log::warn("Reading real-valued array into complex container, consider recompiling with real arithmetic");
    }

    void reset() {
        size_t iline = first_valid_line(m_fname, m_nind, m_indsfirst, m_complex_valued);
        REQUIRE_NE_ALL(iline, ~0ul, "No valid entries found");
        REQUIRE_NE_ALL(m_indsfirst, Tern::Neither, "m_indsfirst is still unspecified");
        REQUIRE_NE_ALL(m_complex_valued, Tern::Neither, "m_complex_valued is still unspecified");
        FileReader::reset(iline);
    }

    virtual bool next(defs::inds &inds, T &v) const {
        std::string line;
        bool result = FileReader::next(line);
        if (!result) return false;
        extract(line, m_nind, m_indsfirst, m_complex_valued, inds, v);
        return true;
    }

    size_t first_valid_line(const std::string &fname, const size_t &nind, Tern &indsfirst, Tern &complex_valued) {
        FileReader iterator(fname);
        std::string line;
        while (iterator.next(line)) {
            if (validate(line, nind, indsfirst, complex_valued)) return iterator.iline();
        }
        return ~0ul; // no valid line found
    }

    static bool validate(const std::string &line, const size_t &nind, Tern &indsfirst, Tern &complex_valued) {
        // can only be complex valued if the real and imag parts are delimited by a comma
        std::vector<double> doubles{};
        auto tokens = string_utils::split(line, "() ,");
        for (const auto& token : tokens) {
            try {
                doubles.push_back(std::stod(token));
            }
            catch (std::invalid_argument &) {
                return false;
            }
        }

        if (complex_valued == Tern::Neither) {
            complex_valued = line.find(',') != ~0ul;
            if (tokens.size() > nind + 1 + complex_valued) return false;
        }
        using namespace float_utils;
        if (indsfirst == Tern::Neither) {
            if (!is_integral(doubles[0]) || (complex_valued && !is_integral(doubles[1]))) {
                // float(s) come first
                indsfirst = false;
                auto begin = doubles.data() + (complex_valued ? 2 : 1);
                auto end = begin + nind;
                return std::all_of(begin, end, is_integral<double>);
            } else if (!is_integral(doubles[doubles.size() - 1]) ||
                       (complex_valued && !is_integral(doubles[doubles.size() - 2]))) {
                // float(s) come last
                indsfirst = true;
                auto begin = doubles.data();
                auto end = begin + nind;
                return std::all_of(begin, end, is_integral<double>);
            }
        }
        return true;
    }

    static void
    extract(const std::string &line, size_t nind, bool indsfirst, bool complex_valued, defs::inds &inds, T &value) {
        const char *p = line.begin().base();
        auto f = [&p, &nind, &complex_valued, &inds, &value](bool indsnext) {
            if (indsnext) {
                for (size_t i = 0ul; i < nind; ++i) {
                    inds[i] = string_utils::read_unsigned(p);
                }
            } else {
                value = string_utils::read_double(p);
                if (consts::is_complex<T>() && complex_valued)
                    complex_utils::set_imag_part(value, string_utils::read_double(p));
            }
        };
        if(indsfirst){ f(1); f(0);}
        else {f(0); f(1);}
    }

    virtual bool inds_in_range(const defs::inds& inds) const {
        return true;
    }

    static bool ind_in_range(size_t ind, size_t extent) {
        return ind < extent || ind==~0ul;
    }
};


#endif //M7_SPARSEARRAYFILEREADER_H
