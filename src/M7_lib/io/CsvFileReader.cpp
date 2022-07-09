//
// Created by Robert J. Anderson on 12/3/21.
//

#include "CsvFileReader.h"
#include "M7_lib/util/String.h"

CsvFileReader::CsvFileReader(const str_t& fname, str_t delimiters, uint_t iline) :
        FileReader(fname, iline), m_delimiters(std::move(delimiters)) {}

bool CsvFileReader::next(strv_t& tokens) {
    if (!FileReader::next(m_work_line)) return false;
    // stop reading at the first empty line:
    if (m_work_line.empty()) return false;
    string::split(m_work_line, tokens, m_delimiters);
    return true;
}

str_t NumericCsvFileReader::c_allowed_chars(".(),+-eE");

bool NumericCsvFileReader::valid_numeric(const str_t& token) {
    for (auto& c: token) {
        if (c >= '0' && c <= '9') continue;
        if (strchr(c_allowed_chars.c_str(), c) != nullptr) continue;
        return false;
    }
    return true;
}

bool NumericCsvFileReader::valid_numeric(const strv_t& tokens) {
    for (auto& token: tokens) if (!valid_numeric(token)) return false;
    return true;
}

NumericCsvFileReader::NumericCsvFileReader(const str_t& fname,
              uint_t ncolumn, str_t delimiters, uint_t iline) :
        CsvFileReader(fname, delimiters, iline), m_ncolumn(ncolumn) {}

bool NumericCsvFileReader::next(strv_t& tokens) {
    do {
        if (!CsvFileReader::next(tokens)) return false;
    } while (tokens.size() != m_ncolumn || !valid_numeric(tokens));
    return true;
}

bool NumericCsvFileReader::parsable_as(const str_t& str, uint_t&) {
    for (auto c : str) if (strchr(c_allowed_chars.c_str(), c) != nullptr) return false;
    return true;
}

bool NumericCsvFileReader::parsable_as(const str_t& str, int&) {
    for (auto c : str) {
        if (c=='-') continue;
        if (strchr(c_allowed_chars.c_str(), c) != nullptr) return false;
    }
    return true;
}

bool NumericCsvFileReader::parsable_as(const str_t&, double&) {
    return true;
}

bool NumericCsvFileReader::parsable_as(const str_t&, float&) {
    return true;
}

void NumericCsvFileReader::parse(const str_t& str, uint_t& v) {
    v = std::stoul(str);
}

void NumericCsvFileReader::parse(const str_t& str, long& v) {
    v = std::stol(str);
}

void NumericCsvFileReader::parse(const str_t& str, int& v) {
    v = std::stoi(str);
}

void NumericCsvFileReader::parse(const str_t& str, double& v) {
    v = std::stod(str);
}

void NumericCsvFileReader::parse(const str_t& str, float& v) {
    v = std::stof(str);
}

uint_t NumericCsvFileReader::ncolumn(const str_t& fname, std::function<uint_t(const str_t&)> iline_fn) {
    if (!FileReader::exists(fname)) return 0ul;
    CsvFileReader file_reader(fname);
    auto iline = iline_fn(fname);
    REQUIRE_NE(iline, ~0ul, "first valid line of CSV data could not be found in file");
    file_reader.skip(iline);
    strv_t tokens;
    file_reader.next(tokens);
    return tokens.size();
}
