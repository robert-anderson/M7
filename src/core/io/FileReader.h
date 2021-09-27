//
// Created by rja on 15/07/2020.
//

#ifndef M7_FILEREADER_H
#define M7_FILEREADER_H

#include <memory>
#include <fstream>

class FileReader {
protected:
    std::unique_ptr<std::ifstream> m_file = nullptr;
    mutable size_t m_iline = ~0ul; // the index of the last extracted line

public:
    const std::string m_fname;
    explicit FileReader(std::string fname, size_t iline = 0ul);

    void reset(size_t iline=0ul);

    virtual ~FileReader();

    const size_t &iline();

    size_t nline();

    bool next(std::string &line) const;

    bool next();

    void skip(size_t nline);

    static std::string to_string(const std::string& fname);

};

#endif //M7_FILEREADER_H
