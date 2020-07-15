//
// Created by rja on 15/07/2020.
//

#ifndef M7_FCIDUMPFILEREADER_H
#define M7_FCIDUMPFILEREADER_H

#include "SparseArrayFileReader.h"
#include <regex>

static std::regex header_terminator_regex{R"(\&END)"};

template<typename T>
class FcidumpFileReader : public SparseArrayFileReader<T> {
    const size_t m_norb;
    const size_t m_isymm;
    const size_t m_nelec;
    const defs::inds m_orbsym;
    const bool m_spin_resolved;
public:
    FcidumpFileReader(const std::string &fname) : SparseArrayFileReader<T>(fname, 4),
    m_norb(read_header_int(fname, "NORB")), m_isymm(m_norb>3?isymm(fname):1),
    m_nelec(read_header_int(fname, "NELEC")),
    m_orbsym(read_header_array(fname, "ORBSYM")),
    m_spin_resolved(read_header_bool(fname, "UHF") || read_header_bool(fname, "TREL"))
            {
    }

    static size_t read_header_int(const std::string &fname, const std::string &label, size_t default_ = 0) {
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
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    static defs::inds read_header_array(const std::string &fname, const std::string &label) {
        return defs::inds{0, 1};
    }

    static size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_ = false) {
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
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }


};

#endif //M7_FCIDUMPFILEREADER_H
