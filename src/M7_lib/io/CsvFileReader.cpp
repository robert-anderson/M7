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

bool NumericCsvFileReader::all_float_parseable(const strv_t& tokens) {
    double tmp;
    for (auto& token: tokens) if (!parse::catching(token, tmp)) return false;
    return true;
}

NumericCsvFileReader::NumericCsvFileReader(const str_t& fname,
              uint_t ncolumn, str_t delimiters, uint_t iline) :
        CsvFileReader(fname, delimiters, iline), m_ncolumn(ncolumn) {}

bool NumericCsvFileReader::next(strv_t& tokens) {
    do {
        if (!CsvFileReader::next(tokens)) return false;
    } while (tokens.size() != m_ncolumn || !all_float_parseable(tokens));
    return true;
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
