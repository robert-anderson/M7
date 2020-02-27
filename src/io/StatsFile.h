//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_STATSFILE_H
#define M7_STATSFILE_H

#include <fstream>
#include <assert.h>
#include <string>
#include <vector>
#include <src/data/Specification.h>
#include <src/parallel/MPIWrapper.h>

struct StatColumn {
    const std::string description;
    const size_t itype;
    std::string output = "none";

    StatColumn(const std::string &description, const size_t itype) : description(description), itype(itype) {}

    template<typename T>
    void write(const T &v) {
        assert(itype == numtypes::itype<T>());
        output = std::to_string(v);
    }

    template<typename T>
    void write(const std::complex<T> &v) {
        assert(itype == numtypes::itype<T>());
        output = std::to_string(v.real()) + " " + std::to_string(v.imag());
    }

    size_t nsubcolumn() const { return numtypes::is_complex(itype) ? 2 : 1; }
};

class StatsFile {
    const std::string m_fname;
    std::ofstream m_file;
    std::vector<StatColumn> m_columns{};

public:
    StatsFile(const std::string &fname="M7.stats") : m_fname(fname), m_file(std::ofstream(fname)) {}

    ~StatsFile() {
        m_file.close();
    }

    template<typename T>
    StatColumn add_column(const std::string &description) {
        assert(numtypes::itype<T>() != ~0ul);
        m_columns.emplace_back(description, numtypes::itype<T>());
        return m_columns.back();
    }

    void write_header();

    void flush();

};


#endif //M7_STATSFILE_H
