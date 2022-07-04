//
// Created by rja on 28/06/22.
//

#include "IntegralReader.h"

void IntegralReader::IterData::update_ranksig_exsig() {
    if (dtype::all_null(m_inds[0], m_inds[1], m_inds[2], m_inds[3])) m_ranksig = 0ul;
    else {
        REQUIRE_FALSE(dtype::any_null(m_inds[0], m_inds[1]),
                      "if not core energy entry, first two inds must be defined");
        m_ranksig = dtype::all_null(m_inds[2], m_inds[3]) ? exsig::ex_single : exsig::ex_double;
    }
    switch (m_ranksig) {
        case 0ul:
            m_exsig = 0ul;
            break;
        case exsig::ex_single:
            m_exsig = m_inds[0] == m_inds[1] ? 0ul : exsig::ex_single;
            break;
        case exsig::ex_double:
            m_exsig = m_inds[0] == m_inds[2] ?
                      (m_inds[1] == m_inds[3] ? 0ul : exsig::ex_single) :
                      (m_inds[1] == m_inds[3] ? exsig::ex_single : exsig::ex_double);
            break;
        default:
            m_exsig = dtype::null(m_exsig);
    }
}

CsvIntegralReader::CsvIntegralReader(const FcidumpInfo& info, bool spin_major) : m_reader(info.m_fname, spin_major){}

bool CsvIntegralReader::next(IntegralReader::IterData& data) {
    auto out = m_reader.next(data.m_inds, data.m_value);
    if (!out) return false;
    data.update_ranksig_exsig();
    if (data.m_ranksig==exsig::ex_double && m_iline_first_2e==~0ul) m_iline_first_2e = m_iline;
    else if (data.m_ranksig==exsig::ex_single && m_iline_first_1e==~0ul) m_iline_first_1e = m_iline;
    ++m_iline;
    return out;
}

void CsvIntegralReader::goto_first_1e() {
    REQUIRE_NE(m_iline_first_1e, ~0ul, "first 1-electron integral yet to be found");
    m_reader.reset(m_iline_first_1e);
    m_iline = m_iline_first_1e;
}

void CsvIntegralReader::goto_first_2e() {
    REQUIRE_NE(m_iline_first_2e, ~0ul, "first 2-electron integral yet to be found");
    m_reader.reset(m_iline_first_2e);
    m_iline = m_iline_first_2e;
}

ham_t CsvIntegralReader::ecore() const {
    return 0.0;
}

bool CsvIntegralReader::complex_valued() const {
    return m_reader.m_complex_valued;
}

Hdf5IntegralReader::Hdf5IntegralReader(const FcidumpInfo& info, Hdf5IntegralReader::KeyNames names, bool) :
        m_reader(info.m_fname), m_names(std::move(names)),
        m_indices_2e(m_reader, m_names.m_2e_inds),
        m_values_2e(m_reader.read_data<v_t<ham_t>>(m_names.m_2e_values)),
        m_indices_1e(m_reader, m_names.m_1e_inds),
        m_values_1e(m_reader.read_data<v_t<ham_t>>(m_names.m_1e_values)){
    REQUIRE_EQ(m_indices_2e.nrow(), m_values_2e.size(),
               "number of 2e matrix index arrays should match the number of values");
    REQUIRE_EQ(m_indices_1e.nrow(), m_values_1e.size(),
               "number of 1e matrix index arrays should match the number of values");
}

bool Hdf5IntegralReader::next(IntegralReader::IterData& data) {
    if (m_iline== m_values_2e.size() + m_values_1e.size()) return false;
    if (m_iline < m_values_2e.size()) {
        const auto ientry = m_iline;
        m_indices_2e.get_row(ientry, data.m_inds);
        integer::shift(data.m_inds, false); // Fortran (1-based) to C (0-based) indexing
        data.update_ranksig_exsig();
        data.m_value = m_values_2e[ientry];
    }
    else {
        const auto ientry = m_iline-m_values_2e.size();
        data.m_inds[0] = m_indices_1e(ientry, 0) - 1;
        data.m_inds[1] = m_indices_1e(ientry, 1) - 1;
        dtype::nullify(data.m_inds[2], data.m_inds[3]);
        data.update_ranksig_exsig();
        data.m_value = m_values_1e[ientry];
    }
    ++m_iline;
    return true;
}

void Hdf5IntegralReader::goto_first_1e() {
    m_iline = m_values_2e.size();
}

void Hdf5IntegralReader::goto_first_2e() {
    m_iline = 0;
}

ham_t Hdf5IntegralReader::ecore() const {
    return m_reader.read_attr<double>(m_names.m_ecore);
}

bool Hdf5IntegralReader::complex_valued() const {
    return false;
}

MolcasHdf5IntegralReader::MolcasHdf5IntegralReader(const FcidumpInfo& info, bool spin_major) :
        Hdf5IntegralReader(info,
                           {"CORE_ENERGY",
                            "FOCK_INDEX",
                            "FOCK_VALUES",
                            "TWO_EL_INT_INDEX",
                            "TWO_EL_INT_VALUES"}, spin_major){}
