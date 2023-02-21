//
// Created by Robert J. Anderson on 12/3/21.
//

#ifndef M7_CSVFILEREADER_H
#define M7_CSVFILEREADER_H

#include <algorithm>

#include <M7_lib/parallel/MPIAssert.h>
#include <M7_lib/util/Parse.h>

#include "FileReader.h"

class CsvFileReader : public FileReader {
protected:
    const str_t m_delimiters;
    str_t m_work_line;
public:

    CsvFileReader(const str_t& fname, str_t delimiters=", ()", uint_t iline=0ul);

    bool next(strv_t& tokens);
};

class NumericCsvFileReader : public CsvFileReader {
    bool all_float_parseable(const strv_t& tokens);

public:
    const uint_t m_ncolumn;
    NumericCsvFileReader(const str_t& fname, uint_t ncolumn, str_t delimiters=", ()", uint_t iline=0ul);

    bool next(strv_t& tokens);

    static uint_t ncolumn(const str_t& fname, std::function<uint_t(const str_t&)> iline_fn);
};

#endif //M7_CSVFILEREADER_H
