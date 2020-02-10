//
// Created by Robert John Anderson on 2020-01-04.
//

#include "FileIterator.h"

FileIterator::FileIterator(const std::string &filename, const size_t &ifirstline):
        m_file(std::ifstream(filename)), m_ifirstline(ifirstline) {
    for (auto i{0ul}; i<m_ifirstline; i++) next();
    assert(m_file.is_open());
}

FileIterator::FileIterator(const std::string &filename):
        FileIterator(filename, 0ul){}

FileIterator::FileIterator(const std::string &filename, const std::regex &regex):
        FileIterator(filename, line_number_from_regex(filename, regex)){}

FileIterator::~FileIterator() {
    m_file.close();
}

const size_t FileIterator::line_number_from_regex(const std::string &filename, const std::regex &regex){
    std::smatch match;
    FileIterator fi(filename);
    std::string line;
    size_t out{0ul};
    while (fi.next(line)) {
        std::regex_search(line, match, regex);
        if (match.size()) return out;
        out++;
    }
    assert(!line.empty());
}

bool FileIterator::next(std::string &line){
    getline(m_file, line);
    return !line.empty();
}

std::string FileIterator::next(){
    std::string line;
    getline(m_file, line);
    return line;
}