//
// Created by rja on 02/04/23.
//

#include "StatsTable.h"

StatsTableBase::StatsTableBase(str_t fname, str_t description, uint_t period, const Row& row) :
        m_fname(std::move(fname)), m_description(std::move(description)), m_period(period),
        m_file(new std::ofstream(m_fname)), m_fields(row.m_fields) {
    write_header();
}

std::map<str_t, uint_t> StatsTableBase::make_format_strings() const {
    std::map<str_t, uint_t> format_strings;
    uint_t nformat = 0ul;
    for (const FieldBase* field: m_fields) {
        auto num_field = dynamic_cast<const NumberFieldBase*>(field);
        auto format_string = num_field->format_string();
        auto it = format_strings.find(format_string);
        if (it == format_strings.end()) format_strings[format_string] = nformat++;
    }
    return format_strings;
}

uint_t StatsTableBase::make_ncolumn() const {
    uint_t ncolumn = 0ul;
    for (const FieldBase* field: m_fields) {
        auto num_field = dynamic_cast<const NumberFieldBase*>(field);
        ncolumn += (num_field->m_is_complex + 1) * num_field->m_nelement;
    }
    return ncolumn;
}


void StatsTableBase::write_header() const {
    if (!m_file) return;
    const auto format_strings = make_format_strings();
    const auto nformat = format_strings.size();
    const auto ncolumn = make_ncolumn();
    *m_file << string::boxed(m_description + " Stats File") <<
            "# Number of statistics output: " << m_fields.size() <<
                 "\n# Number of columns: " << ncolumn <<
                 "\n# Distinct multidimensional formats: " << format_strings.size() << "\n#" <<
                 "\n# Format list (major index first):" << "\n";

    for (uint_t i = 0ul; i < nformat; ++i) {
        for (const auto &it : format_strings) {
            if (it.second == i) {
                *m_file << "# " << static_cast<char>('A' + i) << ": " << it.first << "\n";
                break;
            }
        }
    }
    *m_file << "#\n";

    auto icolumn = 1ul;
    for (const FieldBase *field : m_fields) {
        auto num_field = dynamic_cast<const NumberFieldBase *>(field);
        auto format_string = num_field->format_string();
        auto it = format_strings.find(format_string);
        DEBUG_ASSERT_FALSE(it==format_strings.cend(), "format string should have been found in the map");
        *m_file << "#  " << icolumn << ".  " << num_field->m_name << " (" << static_cast<char>('A' + it->second) << ")\n";
        icolumn += (num_field->m_is_complex + 1) * num_field->m_nelement;
    }
    *m_file << "#\n";
    *m_file << std::flush;
}

str_t StatsTableBase::make_line() const {
    if (!m_file) return "";
    str_t res;
    for (const FieldBase *field : m_fields) {
        auto stats_field = dynamic_cast<const statistic::Base *>(field);
        res += stats_field->stats_string() + " ";
    }
    return res;
}

void StatsTableBase::flush() {
    if (!m_file) return;
    *m_file << make_line() << std::endl;
    for (FieldBase *field : this->m_fields) {
        auto num_field = dynamic_cast<statistic::Base *>(field);
        num_field->reset();
    }
}

void StatsTableBase::commit() {
    if (!m_file) return;
    for (FieldBase *field : this->m_fields)
        dynamic_cast<statistic::Base *>(field)->commit();
    if (!(++m_ncommit % m_period)) flush();
}
