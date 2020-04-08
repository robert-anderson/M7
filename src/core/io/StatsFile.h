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
#include <list>
#include <src/core/table/Table.h>
#include <src/core/table/NumericField.h>
#include <src/core/parallel/MPIWrapper.h>

template<typename T>
class StatsElement : public NumericElement<T> {
    using NumericElement<T>::m_field;
    std::string to_string() const override {
        const T v = NumericElement<T>::to_number();
        return utils::num_to_string(v);
    }
public:

    StatsElement<T>(NumericField<T> *field, char *begin) : NumericElement<T>(field, begin) {}

    StatsElement<T> &operator=(const T &v) override {
        if (!m_field->is_allocated()) m_field->expand_table(1);
        NumericElement<T>::operator=(v);
        return *this;
    }

    StatsElement<T> &operator=(const NumericElement<T> &v) override {
        if (!m_field->is_allocated()) m_field->expand_table(1);
        NumericElement<T>::operator=(v);
        return *this;
    }

    static std::string to_string_fn(const Element *element){
        auto value = *((T *) element->begin());
        if (consts::is_complex<T>())
            return utils::num_to_string(consts::real(value))
                +utils::num_to_string(consts::imag(value));
        else return utils::num_to_string(value);
    }
};

class StatsFile : public Table {
protected:
    const std::string m_fname;
    std::unique_ptr<std::ofstream> m_file;
    size_t m_nflush = 0;
public:
    std::vector<std::function<std::string(const Element*)>> m_printers;

protected:
    void write_header() {
        size_t ncolumn = 0ul;
        for (auto &column : m_fields) ncolumn += column->is_complex();
        *m_file <<
                "################################\n"
                "#    M7 FCIQMC Stats Output    #\n"
                "################################\n"
                "#\n"
                "#\n"
                "# Number of statistics output: " << m_fields.size() <<
                "\n# Number of columns: " << ncolumn << "\n#\n";
        size_t icol = 1ul;
        for (auto &column : m_fields) {
            if (column->is_complex()) {
                *m_file << "#  " << icol++ << ".  " << column->description() << " (real)\n";
                *m_file << "#  " << icol++ << ".  " << column->description() << " (imag)\n";
            } else {
                *m_file << "#  " << icol++ << ".  " << column->description() << " \n";
            }
        }
        *m_file << std::endl;
    }

public:
    StatsFile(const std::string &fname) : Table(), m_fname(fname), m_file(std::make_unique<std::ofstream>(fname)) {
        if (!mpi::i_am_root())
            throw std::runtime_error(
                "StatsFiles must only be instantiated on the root process");
    }

    void flush() {
        if (!m_nflush) {
            if (!is_allocated()) return; //nothing written
            write_header();
        }
        auto printer = m_printers.begin();
        for (auto &column : m_fields) {
            auto element = (*column)(0);
            *m_file << (*printer)(&element) << "  ";
            element.zero();
            printer++;
        }
        *m_file << std::endl;
        m_nflush++;
    }

    ~StatsFile() {
        m_file->close();
    }
};

template<typename T>
class StatsField : public NumericField<T> {
public:
    using NumericField<T>::m_table;
    using NumericField<T>::m_nelement;
    using NumericField<T>::m_element_size;
    using NumericField<T>::begin;

    StatsField(StatsFile *file, size_t nelement = 1, const std::string &description = "") :
        NumericField<T>(file, nelement, description){
        file->m_printers.push_back(StatsElement<T>::to_string_fn);
    }

    StatsElement<T> operator()(const size_t &ielement = 0) {
        assert(ielement < m_nelement);
        if (!m_table->is_allocated()) m_table->expand(1);
        return StatsElement<T>(this, begin(0, 0) + ielement * m_element_size);
    }

    bool is_complex() const override { return consts::is_complex<T>(); }

};

#endif //M7_STATSFILE_H
