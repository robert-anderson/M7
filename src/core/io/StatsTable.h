//
// Created by rja on 12/04/2021.
//

#ifndef M7_STATSTABLE_H
#define M7_STATSTABLE_H

#include "src/core/field/Row.h"
#include "src/core/table/BufferedTable.h"
#include "src/core/table/BufferedFields.h"
#include <fstream>
#include <map>

struct StatsRow : Row {};

namespace statistic {

    struct Base {
        virtual void commit() = 0;

        virtual void reset() = 0;

        virtual std::string stats_string() const = 0;

        template<typename T>
        static std::string stats_string_element(const T &v, size_t denom) {
            return utils::num_to_string(v / denom);
        }

        template<typename T>
        static std::string stats_string_element(const std::complex<T> &v, size_t denom) {
            return utils::num_to_string(v.real() / denom) + " " + utils::num_to_string(v.imag() / denom);
        }
    };

    template<typename T, size_t nind>
    struct Numbers : NdNumberField<T, nind>, Base {
        typedef NdNumberField<T, nind> base_t;
        using base_t::nelement;
        using base_t::operator=;
        using base_t::operator+=;

        const bool m_mean;
        size_t m_ncommit_this_period = 0ul;

        buffered::Numbers<T, nind> m_reduced;

        Numbers(StatsRow *row, NdFormat<nind> format, std::string name = "", bool mean = true) :
                NdNumberField<T, nind>(row, format, name), m_mean(mean), m_reduced(format) {}

        Numbers(const Numbers& other):
                NdNumberField<T, nind>(other), m_mean(other.m_mean), m_reduced(other.m_format){}

        void commit() override {
            if (m_mean) m_reduced += *this;
            else if (!m_ncommit_this_period) m_reduced = *this;
            ++m_ncommit_this_period;
            NumberFieldBase::zero();
        }

        void reset() override {
            m_reduced.zero();
            m_ncommit_this_period = 0ul;
        }

        std::string stats_string() const override {
            std::string tmp;
            for (size_t ielement = 0ul; ielement < nelement(); ++ielement)
                tmp += stats_string_element(m_reduced[ielement], m_mean ? m_ncommit_this_period : 1ul) + " ";
            return tmp;
        }

        Numbers& operator=(const Numbers& other){
            *this = static_cast<const NdNumberField<T, nind>&>(other);
            return *this;
        }
    };

    template<typename T>
    struct Number : Numbers<T, 0ul> {
        typedef Numbers<T, 0ul> base_t;
        using base_t::operator=;
        Number(StatsRow *row, std::string name = "", bool mean=true) : base_t(row, {}, name, mean) {}

        Number(const Number& other): base_t(other){}
        Number& operator=(const Number<T>& other) {return *this;}
    };
}


template<typename row_t>
struct StatsTable : BufferedTable<row_t> {
    static_assert(std::is_base_of<StatsRow, row_t>::value, "Template arg must be derived from StatsRow");
    const std::string m_fname, m_description;
    std::unique_ptr<std::ofstream> m_file;
    const size_t m_period;
    size_t m_ncommit = 0;

    void write_header() const {
        const auto &row = static_cast<const Row &>(Table<row_t>::m_row);
        size_t ncolumn = 0ul;
        std::map<std::string, size_t> m_format_strings;
        size_t nformat = 0ul;
        for (const FieldBase *field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            ncolumn += (num_field->m_is_complex + 1) * num_field->m_nelement;
            auto format_string = num_field->format_string();
            auto it = m_format_strings.find(format_string);
            if (it == m_format_strings.end()) m_format_strings[format_string] = nformat++;
        }
        *m_file << string_utils::boxed(m_description + " Stats File") <<
                "# Number of statistics output: " << row.m_fields.size() <<
                "\n# Number of columns: " << ncolumn <<
                "\n# Distinct multidimensional formats: " << m_format_strings.size() << "\n#" <<
                "\n# Format list (major index first):" << "\n";

        for (size_t i = 0ul; i < nformat; ++i) {
            for (const auto &it : m_format_strings) {
                if (it.second == i) {
                    *m_file << "# " << static_cast<char>('A' + i) << ": " << it.first << "\n";
                    break;
                }
            }
        }
        *m_file << "#\n";

        auto icolumn = 1ul;
        for (const FieldBase *field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            auto format_string = num_field->format_string();
            auto it = m_format_strings.find(format_string);
            if (it == m_format_strings.end()) m_format_strings[format_string] = nformat++;
            *m_file << "#  " << icolumn << ".  " << num_field->m_name <<
                    " (" << static_cast<char>('A' + m_format_strings[num_field->format_string()]) << ")\n";
            icolumn += (num_field->m_is_complex + 1) * num_field->m_nelement;
        }
        *m_file << "#\n";
        *m_file << std::flush;
    }

    StatsTable(std::string fname, std::string description, const row_t &row, size_t period) :
            BufferedTable<row_t>(fname, {row}), m_fname(fname), m_description(description),
            m_file(new std::ofstream(fname)), m_period(period) {
        write_header();
        BufferedTable<row_t>::push_back();
        BufferedTable<row_t>::m_row.restart();
    }

    void commit() {
        const auto &row = static_cast<const Row &>(Table<row_t>::m_row);
        for (FieldBase *field : row.m_fields)
            dynamic_cast<statistic::Base *>(field)->commit();
        if (!(++m_ncommit % m_period)) flush();
    }

private:

    std::string make_line() const {
        std::string res;
        const auto &row = static_cast<const Row &>(Table<row_t>::m_row);
        for (const FieldBase *field : row.m_fields) {
            auto stats_field = dynamic_cast<const statistic::Base *>(field);
            res += stats_field->stats_string() + " ";
        }
        return res;
    }

    void flush() {
        *m_file << make_line() << std::endl;
        auto &row = static_cast<Row &>(Table<row_t>::m_row);
        for (FieldBase *field : row.m_fields) {
            auto num_field = dynamic_cast<statistic::Base *>(field);
            num_field->reset();
        }
    }
};


#endif //M7_STATSTABLE_H
