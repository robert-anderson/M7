//
// Created by rja on 12/04/2021.
//

#ifndef M7_STATSTABLE_H
#define M7_STATSTABLE_H

#include "src/core/field/Row.h"
#include "src/core/table/BufferedTable.h"
#include <fstream>
#include <map>

template<typename row_t>
struct StatsTable : BufferedTable<row_t> {
    const std::string m_fname, m_description;
    std::unique_ptr<std::ofstream> m_file;
    size_t m_nflush = 0;

    void write_header() const {
        const auto& row = static_cast<const Row&>(Table<row_t>::m_row);
        size_t ncolumn = 0ul;
        std::map<std::string, size_t> m_format_strings;
        size_t nformat = 0ul;
        for (const FieldBase* field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            ncolumn+=(num_field->m_is_complex+1)*num_field->m_nelement;
            auto format_string = num_field->format_string();
            auto it = m_format_strings.find(format_string);
            if (it==m_format_strings.end()) m_format_strings[format_string] = nformat++;
        }
        *m_file << string_utils::boxed(m_description+" Stats File") <<
                "# Number of statistics output: " << row.m_fields.size() <<
                "\n# Number of columns: " << ncolumn <<
                "\n# Distinct multidimensional formats: " << m_format_strings.size() << "\n#" <<
                "\n# Format list (major index first):" << "\n";

        for (size_t i=0ul; i<nformat; ++i) {
            for (const auto &it : m_format_strings) {
                if (it.second==i){
                    *m_file << "# " << static_cast<char>('A'+i) << ": " << it.first << "\n";
                    break;
                }
            }
        }
        *m_file << "#\n";

        auto icolumn = 1ul;
        for (const FieldBase* field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            auto format_string = num_field->format_string();
            auto it = m_format_strings.find(format_string);
            if (it==m_format_strings.end()) m_format_strings[format_string] = nformat++;
            *m_file << "#  " << icolumn << ".  " << num_field->m_name <<
                " (" << static_cast<char>('A'+m_format_strings[num_field->format_string()]) << ")\n";
            icolumn+=(num_field->m_is_complex+1)*num_field->m_nelement;
        }
        *m_file << "#\n";
        *m_file << std::flush;
    }

    StatsTable(std::string fname, std::string description, const row_t& row): BufferedTable<row_t>(fname, {row}),
    m_fname(fname), m_description(description), m_file(new std::ofstream(fname)){
        write_header();
        BufferedTable<row_t>::push_back();
        BufferedTable<row_t>::m_row.restart();
    }

    std::string make_line() const {
        std::string res;
        const auto& row = static_cast<const Row&>(Table<row_t>::m_row);
        for (const FieldBase* field : row.m_fields) {
            auto num_field = dynamic_cast<const NumberFieldBase *>(field);
            res+=num_field->stats_string()+" ";
        }
        return res;
    }

    void flush() {
        *m_file << make_line() << std::endl;
        auto& row = static_cast<Row&>(Table<row_t>::m_row);
        for (FieldBase* field : row.m_fields) {
            auto num_field = dynamic_cast<NumberFieldBase *>(field);
            num_field->zero();
        }
        m_nflush++;
    }
};


#endif //M7_STATSTABLE_H
