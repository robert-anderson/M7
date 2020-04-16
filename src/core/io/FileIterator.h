//
// Created by Robert John Anderson on 2020-01-04.
//

#ifndef M7_FILEITERATOR_H
#define M7_FILEITERATOR_H

#include <fstream>
#include <regex>
#include <memory>


class FileIterator {
protected:
    std::unique_ptr<std::ifstream> m_file;
    const size_t m_ifirstline;
public:
    FileIterator(const std::string &filename, const size_t &ifirstline);

    FileIterator(const std::string &filename);

    FileIterator(const std::string &filename, const std::regex &regex);

    virtual ~FileIterator();

    static const size_t line_number_from_regex(const std::string &, const std::regex &);

    bool next(std::string &);

    std::string next();
};


#endif //M7_FILEITERATOR_H
