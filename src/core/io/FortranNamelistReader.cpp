//
// Created by anderson on 12/3/21.
//

#include "FortranNamelistReader.h"

std::regex FortranNamelistReader::c_header_terminator_regex(R"(\&END)");

void FortranNamelistReader::read(size_t &v, const std::string &label, size_t default_) {
    v = default_;
    if (!m_exists) return;
    const auto regex = std::regex(label + R"(\s*\=\s*[0-9]+)");
    FileReader file_reader(m_fname);
    std::string line;
    std::smatch match;
    while (file_reader.next(line)) {
        std::regex_search(line, match, regex);
        if (match.size()) {
            std::string match_string(match.str());
            std::regex_search(match_string, match, std::regex(R"([0-9]+)"));
            v = std::atol(match.str().c_str());
            return;
        }
        std::regex_search(line, match, c_header_terminator_regex);
        if (match.size()) break;
    }
}


void FortranNamelistReader::read(defs::inds &v, const std::string &label, long offset, defs::inds default_) {
    v = default_;
    if (!m_exists) return;
    std::vector<std::string> tokens;
    defs::inds inds;
    const auto regex = std::regex(label + R"(\s*\=\s*[0-9\,\s]+)");
    FileReader file_reader(m_fname);
    std::string line;
    std::smatch match;
    while (file_reader.next(line)) {
        std::regex_search(line, match, regex);
        if (match.size()) {
            // strip away label so any numeric digits in the label are not picked up in error by parse_int_array
            std::string content = match.str();
            std::regex_search(content, match, std::regex(R"(\=.*)"));
            content = match.str();
            string_utils::split(content, tokens, ", ");
            NumericCsvFileReader::parse(tokens.cbegin(), tokens.cend(), inds);
            for (auto &ind: inds) ind += offset;
            return;
        }
        std::regex_search(line, match, c_header_terminator_regex);
        if (match.size()) break;
    }
}

void FortranNamelistReader::read(bool &v, const std::string &label, bool default_) {
    v = default_;
    if (!m_exists) return;
    std::regex regex;
    if (default_) regex = std::regex(label + R"(\s?\=\s?\.FALSE\.)");
    else regex = std::regex(label + R"(\s?\=\s?\.TRUE\.)");
    FileReader file_reader(m_fname);
    std::string line;
    std::smatch match;
    std::string match_string;
    while (file_reader.next(line)) {
        std::regex_search(line, match, regex);
        if (match.size()) {
            v = !default_;
            return;
        }
        std::regex_search(line, match, c_header_terminator_regex);
        if (match.size()) return;
    }
}

size_t FortranNamelistReader::read_int(const std::string &label, size_t default_) {
    size_t tmp;
    read(tmp, label, default_);
    return tmp;
}

defs::inds FortranNamelistReader::read_int_array(const std::string &label, long offset, defs::inds default_) {
    defs::inds tmp;
    read(tmp, label, offset, default_);
    return tmp;
}

bool FortranNamelistReader::read_bool(const std::string &label, size_t default_) {
    bool tmp;
    read(tmp, label, default_);
    return tmp;
}
