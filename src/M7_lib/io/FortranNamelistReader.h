//
// Created by Robert J. Anderson on 12/3/21.
//

#ifndef M7_FORTRANNAMELISTREADER_H
#define M7_FORTRANNAMELISTREADER_H

#include <utility>
#include <regex>

#include "CsvFileReader.h"

/**
 * When dealing with legacy formats such as the FCIDUMP, we need to be able to parse information stored in a header
 * formatted as a Fortran namelist
 */
class FortranNamelistReader {
    const bool m_exists;

public:

    const std::string m_fname;
    explicit FortranNamelistReader(std::string fname);

    static std::string isolate_value(const std::string &line, const std::string &label);

    static bool contains_terminator(const std::string &line);

    static std::vector<std::string> read(const std::string &line, const std::string &label);

    std::vector<std::string> read(const std::string &label) const;

    template<typename T>
    void read(std::vector<T>& v, const std::string &label, long offset=0, std::vector<T> default_ = {}) const {
        v = default_;
        auto tokens = read(label);
        if (tokens.empty()) return;
        NumericCsvFileReader::parse(tokens.cbegin(), tokens.cend(), v);
        for (auto &i: v) i+=offset;
    }

    void read(std::vector<bool>& v, const std::string &label, long offset=0, std::vector<bool> default_ = {}) const;

    std::vector<long> read_ints(const std::string &label, long offset=0, std::vector<long> default_ = {}) const;

    std::vector<size_t> read_uints(const std::string &label, long offset=0, std::vector<size_t> default_ = {}) const;

    std::vector<bool> read_bools(const std::string &label, long offset=0, std::vector<bool> default_ = {}) const;

    long read_int(const std::string &label, long default_ = 0l) const;

    size_t read_uint(const std::string &label, size_t default_ = 0ul) const;

    bool read_bool(const std::string &label, bool default_ = false) const;

};

#endif //M7_FORTRANNAMELISTREADER_H
