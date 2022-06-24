//
// Created by Robert J. Anderson on 27/09/2021.
//

#include <M7_lib/parallel/MPIAssert.h>

#include "FileReader.h"
#include "M7_lib/util/SmartPtr.h"

FileReader::FileReader(std::string fname, uint_t iline) : m_fname(std::move(fname)) {
    reset(iline);
}

void FileReader::reset(uint_t iline) {
    if (m_file) m_file->close();
    m_file = smart_ptr::make_unique<std::ifstream>(m_fname);
    REQUIRE_TRUE(m_file->is_open(), "File not found: " + m_fname);
    m_iline = ~0ul;
    skip(iline);
}

FileReader::~FileReader() {
    m_file->close();
}

const uint_t &FileReader::iline() {
    return m_iline;
}

uint_t FileReader::nline() {
    while (next()){}
    auto nline = iline();
    reset();
    return nline;
}

bool FileReader::next(std::string &line) const {
    m_iline++;
    getline(*m_file, line);
    if (line.empty()) return false;
    // deal with carriage return character
    if (line[line.size()-1]=='\r') line.resize(line.size()-1);
    return !line.empty();
}

bool FileReader::next() {
    std::string tmp;
    return next(tmp);
}

void FileReader::skip(uint_t nline) {
    for (uint_t i = 0ul; i < nline; ++i) next();
}

std::string FileReader::to_string(const std::string &fname) {
    std::string all;
    std::string line;
    FileReader reader(fname);
    while (reader.next(line)) {
        all.append("\n"+line);
        line.clear();
    }
    return all;
}

bool FileReader::exists(const std::string &fname) {
    auto f = smart_ptr::make_unique<std::ifstream>(fname);
    bool res = f->is_open();
    if (res) f->close();
    return res;
}
