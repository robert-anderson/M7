//
// Created by Robert J. Anderson on 12/04/2021.
//

#ifndef M7_STATSTABLE_H
#define M7_STATSTABLE_H

#include <fstream>
#include <map>
#include <utility>

#include <M7_lib/field/Row.h>
#include <M7_lib/table/BufferedTable.h>
#include <M7_lib/table/BufferedFields.h>
#include <M7_lib/util/String.h>

struct StatsRow : Row {};

namespace statistic {

    struct Base {
        static constexpr uint_t c_default_fp = 8ul;
        virtual void commit() = 0;

        virtual void reset() = 0;

        virtual str_t stats_string() const = 0;

        template<typename T>
        static str_t stats_string_element(T v, uint_t denom, uint_t fp=c_default_fp) {
            return convert::to_string(v / denom, fp);
        }

        template<typename T>
        static str_t stats_string_element(std::complex<T> v, uint_t denom, uint_t fp=c_default_fp) {
            return convert::to_string(v.real() / denom, fp) + " " + convert::to_string(v.imag() / denom, fp);
        }
    };

    template<typename T, uint_t nind>
    struct Numbers : NdNumberField<T, nind>, Base {
        typedef NdNumberField<T, nind> base_t;
        using base_t::nelement;
        using base_t::operator=;
        using base_t::operator+=;

        const bool m_mean;
        /**
         * floating point precision
         */
        const uint_t m_fp;
        uint_t m_ncommit_this_period = 0ul;

        v_t<T> m_reduced;

        Numbers(StatsRow *row, NdFormat<nind> format, str_t name = "", bool mean = true, uint_t fp=c_default_fp) :
                NdNumberField<T, nind>(row, format, name), m_mean(mean), m_fp(fp), m_reduced(this->m_nelement) {}

        Numbers(const Numbers& other):
                NdNumberField<T, nind>(other), m_mean(other.m_mean),
                m_fp(other.m_fp), m_reduced(other.m_nelement){}

        void commit() override {
            if (m_mean) this->add_to(m_reduced);
            else if (!m_ncommit_this_period) this->copy_to(m_reduced);
            ++m_ncommit_this_period;
            NumberFieldBase::zero();
        }

        void reset() override {
            m_reduced.assign(this->m_nelement, 0);
            m_ncommit_this_period = 0ul;
        }

        str_t stats_string() const override {
            str_t tmp;
            for (uint_t ielement = 0ul; ielement < nelement(); ++ielement)
                tmp += stats_string_element(m_reduced[ielement], m_mean ? m_ncommit_this_period : 1ul, m_fp) + " ";
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
        Number(StatsRow *row, str_t name = "", bool mean=true) : base_t(row, {}, name, mean) {}

        Number(const Number& other): base_t(other){}
        Number& operator=(const Number<T>& other) {return *this;}
    };
}


template<typename row_t>
struct StatsTable : buffered::Table<row_t> {
    static_assert(std::is_base_of<StatsRow, row_t>::value, "Template arg must be derived from StatsRow");
    const str_t m_fname, m_description;
    std::unique_ptr<std::ofstream> m_file;
    const uint_t m_period;
    uint_t m_ncommit = 0;

    void write_header() const {
        auto& row = this->m_row;
        uint_t ncolumn = 0ul;
        std::map<str_t, uint_t> m_format_strings;
        uint_t nformat = 0ul;
        for (const FieldBase *field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            ncolumn += (num_field->m_is_complex + 1) * num_field->m_nelement;
            auto format_string = num_field->format_string();
            auto it = m_format_strings.find(format_string);
            if (it == m_format_strings.end()) m_format_strings[format_string] = nformat++;
        }
        *m_file << string::boxed(m_description + " Stats File") <<
                "# Number of statistics output: " << row.m_fields.size() <<
                "\n# Number of columns: " << ncolumn <<
                "\n# Distinct multidimensional formats: " << m_format_strings.size() << "\n#" <<
                "\n# Format list (major index first):" << "\n";

        for (uint_t i = 0ul; i < nformat; ++i) {
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

    StatsTable(str_t fname, str_t description, const row_t &row, uint_t period) :
            buffered::Table<row_t>(fname, {row}), m_fname(fname), m_description(std::move(description)),
            m_file(new std::ofstream(fname)), m_period(period) {
        write_header();
        buffered::Table<row_t>::push_back();
        Table<row_t>::m_row.restart();
    }

    void commit() {
        for (FieldBase *field : this->m_row.m_fields)
            dynamic_cast<statistic::Base *>(field)->commit();
        if (!(++m_ncommit % m_period)) flush();
    }

private:

    str_t make_line() const {
        str_t res;
        const auto& row = this->m_row;
        for (const FieldBase *field : row.m_fields) {
            auto stats_field = dynamic_cast<const statistic::Base *>(field);
            res += stats_field->stats_string() + " ";
        }
        return res;
    }

    void flush() {
        *m_file << make_line() << std::endl;
        for (FieldBase *field : this->m_row.m_fields) {
            auto num_field = dynamic_cast<statistic::Base *>(field);
            num_field->reset();
        }
    }
};


#endif //M7_STATSTABLE_H
