//
// Created by Robert J. Anderson on 12/3/21.
//

#include "FortranNamelistReader.h"
#include "M7_lib/util/String.h"

FortranNamelistReader::FortranNamelistReader(str_t fname):
        m_exists(FileReader::exists(fname)), m_fname(std::move(fname)){
    if (!m_exists) REQUIRE_TRUE(m_fname.empty(), "cannot read Fortran namelist header from non-existent file: "+m_fname);
}

str_t FortranNamelistReader::isolate_value(const str_t &line, const str_t &label) {
    auto ibegin = line.find(label + "=");
    if (ibegin >= line.size()) return "";
    // go to end of label and = sign
    ibegin += label.size();
    // strip whitespace before value
    while(line[++ibegin]==' '){}
    auto substring = line.substr(ibegin, str_t::npos);
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

bool FortranNamelistReader::contains_terminator(const str_t &line) {
    if (line.find('/')<line.size()) return true;
    if (line.find("&END")<line.size()) return true;
    if (line.find("$END")<line.size()) return true;
    return false;
}

strv_t FortranNamelistReader::read(const str_t &line, const str_t &label) {
    strv_t tokens;
    auto value = isolate_value(line, label);
    string::split(value, tokens, ", ");
    return tokens;
}

strv_t FortranNamelistReader::read(const str_t &label) const {
    if (!m_exists) return {};
    FileReader file_reader(m_fname);
    str_t line;
    while (file_reader.next(line)) {
        auto tokens = read(line, label);
        if (!tokens.empty()) return tokens;
        if (contains_terminator(line)) break;
    }
    return {};
}

v_t<long> FortranNamelistReader::read_ints(const str_t &label, v_t<long> default_) const {
    v_t<long> tmp;
    read(tmp, label, default_);
    return tmp;
}

v_t<uint_t>
FortranNamelistReader::read_uints(const str_t &label, v_t<uint_t> default_) const {
    v_t<uint_t> tmp;
    read(tmp, label, default_);
    return tmp;
}

v_t<bool> FortranNamelistReader::read_bools(const str_t &label, v_t<bool> default_) const {
    v_t<bool> tmp;
    read(tmp, label, default_);
    return tmp;
}

long FortranNamelistReader::read_int(const str_t &label, long default_) const {
    auto tmp = read_ints(label, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

uint_t FortranNamelistReader::read_uint(const str_t &label, uint_t default_) const {
    auto tmp = read_uints(label, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

bool FortranNamelistReader::read_bool(const str_t &label, bool default_) const {
    auto tmp = read_bools(label, {default_});
    REQUIRE_EQ(tmp.size(), 1ul, "found more than one value for scalar");
    return tmp[0];
}

void FortranNamelistReader::read(v_t<bool> &v, const str_t &label, v_t<bool> default_) const {
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
