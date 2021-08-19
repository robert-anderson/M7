//
// Created by rja on 19/08/2021.
//

#ifndef M7_HAMILTONIANFILEREADER_H
#define M7_HAMILTONIANFILEREADER_H

#include <regex>
#include "SparseArrayFileReader.h"

/**
 * base class for all Hamiltonian-defining files with a FORTRAN namelist header, and all values defined at the beginning
 * of the line with indices following
 */
struct HamiltonianFileReader : public SparseArrayFileReader<defs::ham_t> {
    typedef SparseArrayFileReader<defs::ham_t> base_t;
    const std::regex m_header_terminator_regex;
    const size_t m_norb;
    const bool m_spin_resolved;
    const size_t m_nspatorb;

    HamiltonianFileReader(const std::string &fname, size_t nind, bool spin_major):
    base_t(fname, nind, false),
    m_header_terminator_regex(R"(\&END)"), m_norb(read_header_int(fname, "NORB")),
    m_spin_resolved(read_header_bool(fname, "UHF") || read_header_bool(fname, "TREL")),
    m_nspatorb(m_spin_resolved?m_norb/2:m_norb){}

    size_t read_header_int(const std::string &fname, const std::string &label, size_t default_=0) {
        const auto regex = std::regex(label + R"(\s*\=\s*[0-9]+)");
        FileReader iterator(fname);
        std::string line;
        std::smatch match;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) {
                std::string match_string(match.str());
                std::regex_search(match_string, match, std::regex(R"([0-9]+)"));
                return std::atol(match.str().c_str());
            }
            std::regex_search(line, match, m_header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    defs::inds read_header_array(const std::string &fname, const std::string &label) {
        return defs::inds{0, 1};
    }

    size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_=false) {
        std::regex regex;
        if (default_) regex = std::regex(label + R"(\s?\=\s?\.FALSE\.)");
        else regex = std::regex(label + R"(\s?\=\s?\.TRUE\.)");
        FileReader iterator(fname);
        std::string line;
        std::smatch match;
        std::string match_string;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) return !default_;
            std::regex_search(line, match, m_header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

};


#endif //M7_HAMILTONIANFILEREADER_H
