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
    StatsFile(const std::string &fname) : m_fname(fname), m_file(std::ofstream(fname)) {}

    ~StatsFile() {
        m_file.close();
    }

    template<typename T>
    StatColumn add_column(const std::string &description) {
        assert(numtypes::itype<T>() != ~0ul);
        m_columns.emplace_back(description, numtypes::itype<T>());
        return m_columns.back();
    }

    void write_header() {
        size_t ncolumn = 0ul;
        for (auto &column : m_columns) ncolumn += column.nsubcolumn();
        m_file <<
               "################################\n"
               "#    M7 FCIQMC Stats Output    #\n"
               "################################\n"
               "#\n"
               "# Number of statistics output: " << m_columns.size() << "\n"
               "# Number of columns: " << ncolumn << "\n#\n";
        size_t icol = 1ul;
        for (auto &column : m_columns) {
            if (numtypes::is_complex(column.itype)) {
                m_file << "#  " << icol++ << ".  " << column.description << " (real)\n";
                m_file << "#  " << icol++ << ".  " << column.description << " (imag)\n";
            } else {
                m_file << "#  " << icol++ << ".  " << column.description << " \n";
            }
        }
        m_file << std::endl;
    }

    void flush() {
        for (auto &column : m_columns) {
            m_file << column.output << "  ";
            column.output = "none";
        }
        m_file << std::endl;
    }

};


#endif //M7_STATSFILE_H
