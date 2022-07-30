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

    const str_t m_fname;
    explicit FortranNamelistReader(str_t fname);

    static str_t isolate_value(const str_t &line, const str_t &label);

    static bool contains_terminator(const str_t &line);

    static strv_t read(const str_t &line, const str_t &label);

    strv_t read(const str_t &label) const;

    template<typename T>
    void read(v_t<T>& v, const str_t &label, v_t<T> default_ = {}) const {
        v = default_;
        auto tokens = read(label);
        if (tokens.empty()) return;
        NumericCsvFileReader::parse(tokens.cbegin(), tokens.cend(), v);
    }

    void read(v_t<bool>& v, const str_t &label, v_t<bool> default_ = {}) const;

    v_t<long> read_ints(const str_t &label, v_t<long> default_ = {}) const;

    v_t<uint_t> read_uints(const str_t &label, v_t<uint_t> default_ = {}) const;

    v_t<bool> read_bools(const str_t &label, v_t<bool> default_ = {}) const;

    long read_int(const str_t &label, long default_ = 0l) const;

    uint_t read_uint(const str_t &label, uint_t default_ = 0ul) const;

    bool read_bool(const str_t &label, bool default_ = false) const;

};

#endif //M7_FORTRANNAMELISTREADER_H
