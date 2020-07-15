//
// Created by rja on 15/07/2020.
//

#ifndef M7_FILEREADER_H
#define M7_FILEREADER_H

#include <memory>
#include <fstream>

class FileReader {
protected:
    std::unique_ptr<std::ifstream> m_file;
    size_t m_iline = ~0ul; // the index of the last extracted line
public:
    FileReader(const std::string &fname, size_t iline = 0ul) : m_file(new std::ifstream(fname)) {
        if (!m_file->is_open()) throw std::runtime_error("File \"" + fname + "\" not found.");
        skip(iline);
    }

    virtual ~FileReader() {
        m_file->close();
    }

    const size_t &iline() {
        return m_iline;
    }

    bool next(std::string &line) {
        m_iline++;
        getline(*m_file, line);
        return !line.empty();
    }

    bool next() {
        std::string tmp;
        return next(tmp);
    }

    void skip(size_t niter) {
        for (size_t i = 0ul; i < niter; ++i) next();
    }
};

#endif //M7_FILEREADER_H
