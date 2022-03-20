//
// Created by anderson on 12/3/21.
//

#ifndef M7_FORTRANNAMELISTREADER_H
#define M7_FORTRANNAMELISTREADER_H

#include <utility>
#include <regex>

#include <M7_lib/util/utils.h>

#include "CsvFileReader.h"

/**
 * When dealing with legacy formats such as the FCIDUMP, we need to be able to parse information stored in a header
 * formatted as a Fortran namelist
 */
class FortranNamelistReader {
    const bool m_exists;

public:
    static const std::string c_header_terminator;
    const std::string m_fname;
    FortranNamelistReader(std::string fname);

    static std::string isolate_value(const std::string &line, const std::string &label);

    static std::vector<std::string> read(const std::string &line, const std::string &label);

    std::vector<std::string> read(const std::string &label);

    template<typename T>
    void read(std::vector<T>& v, const std::string &label, long offset=0, std::vector<T> default_={}) {
        v = default_;
        auto tokens = read(label);
        if (tokens.empty()) return;
        NumericCsvFileReader::parse(tokens.cbegin(), tokens.cend(), v);
        for (auto &i: v) i+=offset;
    }

    void read(std::vector<bool>& v, const std::string &label, long offset=0, std::vector<bool> default_={});

    std::vector<long> read_ints(const std::string &label, long offset=0, std::vector<long> default_={});

    std::vector<size_t> read_uints(const std::string &label, long offset=0, std::vector<size_t> default_={});

    std::vector<bool> read_bools(const std::string &label, long offset=0, std::vector<bool> default_={});

    long read_int(const std::string &label, long default_=0l);

    size_t read_uint(const std::string &label, size_t default_=0l);

    bool read_bool(const std::string &label, bool default_=0l);

};

#endif //M7_FORTRANNAMELISTREADER_H
