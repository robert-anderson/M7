//
// Created by rja on 15/07/2020.
//

#ifndef M7_FILEREADER_H
#define M7_FILEREADER_H

#include <memory>
#include <fstream>

class FileReader {
protected:
    const std::string m_fname;
    std::unique_ptr<std::ifstream> m_file = nullptr;
    mutable size_t m_iline = ~0ul; // the index of the last extracted line
public:
    FileReader(const std::string &fname, size_t iline = 0ul) : m_fname(fname) {
        std::cout << "File \"" << fname << "\" opened for reading" << std::endl;
        reset(iline);
    }

    void reset(size_t iline=0ul){
        if (m_file) m_file->close();
        m_file = std::unique_ptr<std::ifstream>(new std::ifstream(m_fname));
        if (!m_file->is_open()) throw std::runtime_error("File \"" + m_fname + "\" not found.");
        m_iline = ~0ul;
        skip(iline);
    }

    virtual ~FileReader() {
        m_file->close();
    }

    const size_t &iline() {
        return m_iline;
    }

    bool next(std::string &line) const {
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
