//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_STATSFILE_H
#define M7_STATSFILE_H

#include <fstream>
#include <assert.h>
#include <string>
#include <memory>
#include <utility>
#include <vector>
#include <src/data/Specification.h>
#include <src/parallel/MPIWrapper.h>
#include <list>

struct StatColumn {
protected:
    std::string m_description;
    size_t m_itype;
    std::string m_output;
public:
    explicit StatColumn(std::string description="inactive", const size_t itype=~0ul) :
    m_description(std::move(description)), m_itype(itype) {
        reset();
    }

    template<typename T>
    void write(const T &v) {
        if (m_itype==~0ul)
            throw std::runtime_error("Cannot write output to an inactive StatColumn.");
        if (m_itype != numtypes::itype<T>())
            throw std::runtime_error("Passed type does not match StatColumn type.");
        m_output = std::to_string(v);
    }

    template<typename T>
    void write(const std::complex<T> &v) {
        if (m_itype==~0ul)
            throw std::runtime_error("Cannot write output to an inactive StatColumn.");
        if (m_itype!=numtypes::itype<std::complex<T>>())
            throw std::runtime_error("Passed type does not match StatColumn type.");
        m_output = std::to_string(v.real()) + " " + std::to_string(v.imag());
    }

    size_t nsubcolumn() const { return numtypes::is_complex(m_itype) ? 2 : 1; }

    const std::string &description() const {
        return m_description;
    }

    size_t itype() const {
        return m_itype;
    }

    const std::string &output() const {
        return m_output;
    }

    void reset() {
        m_output = "none";
    }
};

struct StatsFile {
protected:
    const std::string m_fname;
    std::unique_ptr<std::ofstream> m_file;
    std::list<StatColumn> m_columns{};

    StatsFile(const std::string &fname) : m_fname(fname), m_file(std::make_unique<std::ofstream>(fname)) {
        if (!mpi::i_am_root())
            throw std::runtime_error("StatsFiles must only be instantiated on the root process");
    }

public:
    ~StatsFile() {
        m_file->close();
    }

    template<typename T>
    StatColumn* add_column(const std::string &description) {
        assert(numtypes::itype<T>() != ~0ul);
        m_columns.emplace_back(description, numtypes::itype<T>());
        return &m_columns.back();
    }

    void write_header();

    void flush();

};


#endif //M7_STATSFILE_H
