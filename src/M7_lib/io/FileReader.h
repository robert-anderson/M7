//
// Created by Robert J. Anderson on 15/07/2020.
//

#ifndef M7_FILEREADER_H
#define M7_FILEREADER_H

#include <memory>
#include <fstream>
#include "M7_lib/defs.h"

using namespace defs;

class FileReader {
protected:
    std::unique_ptr<std::ifstream> m_file = nullptr;
    mutable uint_t m_iline = ~0ul; // the index of the last extracted line

public:
    const std::string m_fname;
    explicit FileReader(std::string fname, uint_t iline = 0ul);

    void reset(uint_t iline=0ul);

    virtual ~FileReader();

    const uint_t &iline();

    uint_t nline();

    bool next(std::string &line) const;

    bool next();

    void skip(uint_t nline);

    static std::string to_string(const std::string& fname);

    static bool exists(const std::string& fname);
};

#endif //M7_FILEREADER_H
