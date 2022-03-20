//
// Created by anderson on 12/3/21.
//

#include "FortranNamelistReader.h"

const std::string FortranNamelistReader::c_header_terminator = "&END";

FortranNamelistReader::FortranNamelistReader(std::string fname):
        m_exists(FileReader::exists(fname)), m_fname(std::move(fname)){
    if (!m_exists)
        REQUIRE_TRUE(m_fname.empty(), "cannot read Fortran namelist header from non-existent file: "+m_fname);
}

std::string FortranNamelistReader::isolate_value(const std::string &line, const std::string &label) {
    auto ibegin = line.find(label + "=");
    if (ibegin >= line.size()) return "";
    // go to end of label and = sign
    ibegin += label.size();
    // strip whitespace before value
    while(line[++ibegin]==' '){}
    auto substring = line.substr(ibegin, std::string::npos);
    auto iend = substring.find('=');
    if (iend > substring.size()) {
        iend = substring.size();
    }
    else {
        // another key value pair was found after this one
        while(--iend) {
            // walk back to last occurrence of the delimiter before the '='
            if (substring[iend]==',') break;
        }
    }
    return substring.substr(0, iend);
}


std::vector<std::string> FortranNamelistReader::read(const std::string &line, const std::string &label) {
    std::vector<std::string> tokens;
    auto value = isolate_value(line, label);
    string_utils::split(value, tokens, ", ");
    return tokens;
}

std::vector<std::string> FortranNamelistReader::read(const std::string &label) {
    if (!m_exists) return {};
    FileReader file_reader(m_fname);
    std::string line;
    while (file_reader.next(line)) {
        auto tokens = read(line, label);
        if (!tokens.empty()) return tokens;
        if (line.find(c_header_terminator)<line.size()) break;
    }
    return {};
}

std::vector<long> FortranNamelistReader::read_ints(const std::string &label, long offset, std::vector<long> default_) {
    std::vector<long> tmp;
    read(tmp, label, offset, default_);
    return tmp;
}

std::vector<size_t>
FortranNamelistReader::read_uints(const std::string &label, long offset, std::vector<size_t> default_) {
    std::vector<size_t> tmp;
    read(tmp, label, offset, default_);
    return tmp;
}

std::vector<bool> FortranNamelistReader::read_bools(const std::string &label, long offset, std::vector<bool> default_) {
    std::vector<bool> tmp;
    read(tmp, label, offset, default_);
    return tmp;
}

long FortranNamelistReader::read_int(const std::string &label, long default_) {
    auto tmp = read_ints(label, 0, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

size_t FortranNamelistReader::read_uint(const std::string &label, size_t default_) {
    auto tmp = read_uints(label, 0, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

bool FortranNamelistReader::read_bool(const std::string &label, bool default_) {
    auto tmp = read_bools(label, 0, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

void FortranNamelistReader::read(std::vector<bool> &v, const std::string &label,
                                 long offset, std::vector<bool> default_) {
    v = default_;
    auto tokens = read(label);
    if (tokens.empty()) return;
    v.clear();
    for (auto &token: tokens) {
        if (token==".TRUE.") v.push_back(true);
        else if (token==".FALSE.") v.push_back(false);
        else {
            v = default_;
            return;
        }
    }
}
