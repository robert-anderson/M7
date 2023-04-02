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
        static constexpr bool default_float_scientific() {return true;}
        static constexpr uint_t default_float_precision() {return 8ul;}
        static constexpr convert::FloatFmt default_fmt() {
            return {default_float_scientific(), default_float_precision()};
        };
        virtual void commit() = 0;

        virtual void reset() = 0;

        virtual str_t stats_string() const = 0;

        template<typename T>
        static str_t stats_string_element(T v, uint_t denom, convert::FloatFmt fmt=default_fmt()) {
            return convert::to_string(v / T(denom), fmt);
        }

        template<typename T>
        static str_t stats_string_element(std::complex<T> v, uint_t denom, convert::FloatFmt fmt=default_fmt()) {
            const arith::comp_t<T> d = denom;
            return convert::to_string(v.real() / d, fmt) + " " + convert::to_string(v.imag() / d, fmt);
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
         * floating point format
         */
        const convert::FloatFmt m_fmt;
        uint_t m_ncommit_this_period = 0ul;

        v_t<T> m_reduced;

        Numbers(StatsRow *row, NdFormat<nind> format, str_t name = "",
                bool mean = true, convert::FloatFmt fmt=default_fmt()) :
                NdNumberField<T, nind>(row, format, name),
                m_mean(mean), m_fmt(fmt), m_reduced(this->m_nelement) {}

        Numbers(const Numbers& other):
                NdNumberField<T, nind>(other), m_mean(other.m_mean),
                m_fmt(other.m_fmt), m_reduced(other.m_nelement){}

        void commit() override {
            if (m_mean) this->add_to(m_reduced);
            else if (!m_ncommit_this_period) this->copy_to(m_reduced);
            ++m_ncommit_this_period;
            NumberFieldBase::clear();
        }

        void reset() override {
            m_reduced.assign(this->m_nelement, 0);
            m_ncommit_this_period = 0ul;
        }

        str_t stats_string() const override {
            str_t tmp;
            for (uint_t ielement = 0ul; ielement < nelement(); ++ielement)
                tmp += stats_string_element(m_reduced[ielement], m_mean ? m_ncommit_this_period : 1ul, m_fmt) + " ";
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
        Number(StatsRow *row, str_t name = "", bool mean = true, convert::FloatFmt fmt=base_t::default_fmt()) :
            base_t(row, {}, name, mean, fmt) {}

        Number(const Number& other): base_t(other){}
        Number& operator=(const Number<T>& other) {return *this;}
    };
}

struct StatsTableBase {
    const str_t m_fname, m_description;
    const uint_t m_period;
    std::unique_ptr<std::ofstream> m_file;
    const v_t<FieldBase *> m_fields;
    uint_t m_ncommit = 0;

    StatsTableBase(str_t fname, str_t description, uint_t period, const Row& row);

private:

    std::map<str_t, uint_t> make_format_strings() const;

    uint_t make_ncolumn() const;

    void write_header() const;

    str_t make_line() const;

    void flush();

public:
    void commit() ;

};

template<typename row_t>
struct StatsTable : public buffered::Table<row_t>, StatsTableBase {
    static_assert(std::is_base_of<StatsRow, row_t>::value, "Template arg must be derived from StatsRow");

    StatsTable(str_t fname, str_t description, const row_t &row, uint_t period) :
            buffered::Table<row_t>(fname, {row}),
            StatsTableBase(fname, std::move(description), period, this->m_row){
        buffered::Table<row_t>::push_back();
        Table<row_t>::m_row.restart();
    }
};


#endif //M7_STATSTABLE_H
