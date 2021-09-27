//
// Created by rja on 27/09/2021.
//

#include <src/core/parallel/MPIAssert.h>
#include "FileReader.h"

FileReader::FileReader(std::string fname, size_t iline) : m_fname(std::move(fname)) {
    reset(iline);
}

void FileReader::reset(size_t iline) {
    if (m_file) m_file->close();
    m_file = std::unique_ptr<std::ifstream>(new std::ifstream(m_fname));
    REQUIRE_TRUE(m_file->is_open(), "File not found: {}" + m_fname);
    m_iline = ~0ul;
    skip(iline);
}

FileReader::~FileReader() {
    m_file->close();
}

const size_t &FileReader::iline() {
    return m_iline;
}

size_t FileReader::nline() {
    while (next()){}
    auto nline = iline();
    reset();
    return nline;
}

bool FileReader::next(std::string &line) const {
    m_iline++;
    getline(*m_file, line);
    return !line.empty();
}

bool FileReader::next() {
    std::string tmp;
    return next(tmp);
}

void FileReader::skip(size_t nline) {
    for (size_t i = 0ul; i < nline; ++i) next();
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