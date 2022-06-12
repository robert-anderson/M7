//
// Created by Robert J. Anderson on 12/3/21.
//

#include "CsvFileReader.h"
#include "M7_lib/util/String.h"

CsvFileReader::CsvFileReader(const std::string &fname, std::string delimiters, size_t iline) :
        FileReader(fname, iline), m_delimiters(std::move(delimiters)) {}

bool CsvFileReader::next(std::vector<std::string> &tokens) {
    if (!FileReader::next(m_work_line)) return false;
    utils::string::split(m_work_line, tokens, m_delimiters);
    return true;
}

std::string NumericCsvFileReader::c_allowed_chars(".(),+-eE");

bool NumericCsvFileReader::valid_numeric(const std::string &token) {
    for (auto &c: token) {
        if (c >= '0' && c <= '9') continue;
        if (strchr(c_allowed_chars.c_str(), c) != nullptr) continue;
        return false;
    }
    return true;
}

bool NumericCsvFileReader::valid_numeric(const std::vector<std::string> &tokens) {
    for (auto &token: tokens) if (!valid_numeric(token)) return false;
    return true;
}

NumericCsvFileReader::NumericCsvFileReader(const std::string &fname,
              size_t ncolumn, std::string delimiters, size_t iline) :
        CsvFileReader(fname, delimiters, iline), m_ncolumn(ncolumn) {}

bool NumericCsvFileReader::next(std::vector<std::string> &tokens) {
    do {
        if (!CsvFileReader::next(tokens)) return false;
    } while (tokens.size() != m_ncolumn || !valid_numeric(tokens));
    return true;
}

bool NumericCsvFileReader::parsable_as(const std::string &str, size_t &v) {
    for (auto c : str) if (strchr(c_allowed_chars.c_str(), c) != nullptr) return false;
    return true;
}

bool NumericCsvFileReader::parsable_as(const std::string &str, int &v) {
    for (auto c : str) {
        if (c=='-') continue;
        if (strchr(c_allowed_chars.c_str(), c) != nullptr) return false;
    }
    return true;
}

bool NumericCsvFileReader::parsable_as(const std::string &str, double &v) {
    return true;
}

bool NumericCsvFileReader::parsable_as(const std::string &str, float &v) {
    return true;
}

void NumericCsvFileReader::parse(const std::string &str, size_t &v) {
    v = std::stoul(str);
}

void NumericCsvFileReader::parse(const std::string &str, long &v) {
    v = std::stol(str);
}

void NumericCsvFileReader::parse(const std::string &str, int &v) {
    v = std::stoi(str);
}

void NumericCsvFileReader::parse(const std::string &str, double &v) {
    v = std::stod(str);
}

void NumericCsvFileReader::parse(const std::string &str, float &v) {
    v = std::stof(str);
}

size_t NumericCsvFileReader::ncolumn(const std::string &fname, std::function<size_t(const std::string &)> iline_fn) {
    if (!FileReader::exists(fname)) return 0ul;
    CsvFileReader file_reader(fname);
    auto iline = iline_fn(fname);
    REQUIRE_NE(iline, ~0ul, "first valid line of CSV data could not be found in file");
    file_reader.skip(iline);
    std::vector<std::string> tokens;
    file_reader.next(tokens);
    return tokens.size();
}
