//
// Created by anderson on 12/3/21.
//

#ifndef M7_FORTRANNAMELISTREADER_H
#define M7_FORTRANNAMELISTREADER_H

#include <utility>
#include <regex>
#include <src/core/util/utils.h>
#include "CsvFileReader.h"

/**
 * When dealing with legacy formats such as the FCIDUMP, we need to be able to parse information stored in a header
 * formatted as a Fortran namelist
 */
class FortranNamelistReader {
    const bool m_exists;

public:
    static std::regex c_header_terminator_regex;
    const std::string m_fname;
    FortranNamelistReader(std::string fname):
        m_exists(FileReader::exists(fname)), m_fname(std::move(fname)){}

    void read(size_t& v, const std::string &label, size_t default_=0);

    void read(defs::inds& v, const std::string &label, long offset=0, defs::inds default_={});

    void read(bool& v, const std::string &label, bool default_=false);

    size_t read_int(const std::string &label, size_t default_=0);

    defs::inds read_int_array(const std::string &label, long offset=0, defs::inds default_={});

    bool read_bool(const std::string &label, size_t default_=false);
};

#endif //M7_FORTRANNAMELISTREADER_H
