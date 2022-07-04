//
// Created by Robert J. Anderson on 15/07/2020.
//

#ifndef M7_FILEREADER_H
#define M7_FILEREADER_H

#include <memory>
#include <fstream>
#include "M7_lib/defs.h"


class FileReader {
protected:
    std::unique_ptr<std::ifstream> m_file = nullptr;
    mutable uint_t m_iline = ~0ul; // the index of the last extracted line

public:
    const str_t m_fname;
    explicit FileReader(str_t fname, uint_t iline = 0ul);

    void reset(uint_t iline=0ul);

    virtual ~FileReader();

    const uint_t &iline();

    uint_t nline();

    bool next(str_t &line) const;

    bool next();

    void skip(uint_t nline);

    static str_t to_string(const str_t& fname);

    static bool exists(const str_t& fname);
};

#endif //M7_FILEREADER_H
