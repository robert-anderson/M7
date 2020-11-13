//
// Created by Robert John Anderson on 2020-02-23.
//

#ifndef M7_STATSFILE_H
#define M7_STATSFILE_H

#include <fstream>
#include <memory>
#include "StatsColumn.h"
#include "src/core/util/utils.h"

struct StatsSpecifier {
    std::string m_description;
    std::vector<StatsColumnBase *> m_columns;

    StatsSpecifier(std::string description): m_description(description){}
    void add_column(StatsColumnBase *column) {
        m_columns.push_back(column);
    }
};


template<typename spec_t>
struct StatsFile : spec_t {
    static_assert(std::is_base_of<StatsSpecifier, spec_t>::value, "Template arg must be derived from StatsSpecifier");
    const std::string m_fname;
    std::unique_ptr<std::ofstream> m_file;
    size_t m_nflush = 0;
    typedef std::unique_ptr<StatsFile<spec_t>> ptr_t;

    template<typename ...Args>
    StatsFile(std::string fname, Args... spec_args):
            spec_t(spec_args...), m_fname(fname), m_file(new std::ofstream(fname)) {
                write_header();
            }

    using spec_t::m_columns;

    const std::vector<StatsColumnBase *> &columns() const {
        return m_columns;
    }

    void write_header() const {
        size_t ncolumn = 0ul;
        for (const StatsColumnBase * column : columns())
            ncolumn += column->m_nelement*column->m_nsubcolumn;
        auto description = static_cast<const StatsSpecifier*>(this)->m_description;
        *m_file << string_utils::boxed(description+" Stats File") <<
                "# Number of statistics output: " << columns().size() <<
                "\n# Number of columns: " << ncolumn << "\n#\n";
        size_t icol = 1ul;
        for (const StatsColumnBase * column : columns()) {
            auto format_enum = column->format_enum();
            defs::inds inds(column->m_shape.size());
            while (format_enum.next(inds)) {
                auto shape_string = inds.empty() ? "" : utils::to_string(inds);
                if (column->m_nsubcolumn == 2) {
                    *m_file << "#  " << icol++ << ".  " << column->m_description <<
                            " " << shape_string << " (real)\n";
                    *m_file << "#  " << icol++ << ".  " << column->m_description <<
                            " " << shape_string << " (imag)\n";
                } else {
                    *m_file << "#  " << icol++ << ".  " << column->m_description <<
                            " " << shape_string<< "\n";
                }
            }
        }
        *m_file << std::flush;
    }

    std::string to_string() const {
        std::string res;
        for (const StatsColumnBase * column : columns()) res+=column->to_string();
        return res;
    }

    void flush() {
        *m_file << to_string() << std::endl;
        for (StatsColumnBase *column : m_columns) column->zero();
        m_nflush++;
    }
};

#if 0
    std::pair<T,T> mean_std(size_t istart, size_t iend) const {
        ASSERT(iend>istart);
        auto cbegin = m_series.cbegin(); std::advance(cbegin, istart);
        auto cend = m_series.cbegin(); std::advance(cend, iend);
        return stat_utils::mean_std<T>(cbegin, cend);
    }

    std::pair<T,T> mean_std(size_t istart) {
        return mean_std(istart, m_series.size());
    }

    std::pair<T,T> mean_std() {
        return mean_std(0);
    }

#endif //M7_STATSFILE_H
#endif //M7_STATSFILE_H
