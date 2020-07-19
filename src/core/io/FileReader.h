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
    const size_t m_nline;
    mutable size_t m_iline = ~0ul; // the index of the last extracted line

    size_t calc_nline(){
        while (next()){}
        auto nline = iline();
        reset();
        return nline;
    }
public:
    FileReader(const std::string &fname, size_t iline = 0ul) : m_fname(fname), m_nline(calc_nline()) {
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
    const size_t &nline() const {
        return m_nline;
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

    void skip(size_t nline) {
        for (size_t i = 0ul; i < nline; ++i) next();
    }

};

#endif //M7_FILEREADER_H
