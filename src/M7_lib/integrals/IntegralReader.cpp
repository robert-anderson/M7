//
// Created by rja on 28/06/22.
//

#include "IntegralReader.h"

void IntegralReader::IterData::update_ranksig_exsig() {
    if (dtype::all_null(m_inds[0], m_inds[1], m_inds[2], m_inds[3])) m_ranksig = opsig::c_zero;
    else if (!dtype::is_null(m_inds[0]) && dtype::is_null(m_inds[1])) {
        // this is an orbital energy, not a Hamiltonian coefficient, so discard
        m_ranksig = opsig::c_invalid;
        m_exsig = opsig::c_invalid;
        return;
    }
    else {
        REQUIRE_FALSE(dtype::is_null(m_inds[0]),
                      "if not core energy entry, at least the first index must be validly defined");
        m_ranksig = dtype::all_null(m_inds[2], m_inds[3]) ? opsig::c_sing : opsig::c_doub;
    }
    switch (m_ranksig) {
        case opsig::c_zero:
            m_exsig = opsig::c_zero;
            break;
        case opsig::c_sing:
            m_exsig = m_inds[0] == m_inds[1] ? opsig::c_zero : opsig::c_sing;
            break;
        case opsig::c_doub:
            m_exsig = m_inds[0] == m_inds[2] ?
                      (m_inds[1] == m_inds[3] ? opsig::c_zero : opsig::c_sing) :
                      (m_inds[1] == m_inds[3] ? opsig::c_sing : opsig::c_doub);
            break;
        default:
            m_exsig = opsig::c_invalid;
    }
}

CsvIntegralReader::CsvIntegralReader(const FcidumpInfo& info) : m_reader(info){}

bool CsvIntegralReader::next(IntegralReader::IterData& data) {
    auto out = m_reader.next(data.m_inds, data.m_value);
    if (!out) return false;
    data.update_ranksig_exsig();
    if ((data.m_ranksig==opsig::c_doub) && m_iline_first_2e==~0ul) m_iline_first_2e = m_iline;
    else if ((data.m_ranksig==opsig::c_sing) && m_iline_first_1e==~0ul) m_iline_first_1e = m_iline;
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

bool CsvIntegralReader::spin_conserving(uint_t iex) const {
    REQUIRE_TRUE((iex==1ul) || (iex==2ul), "integrals only refer to 1 or 2 body Hamiltonian coefficients");
    switch (iex) {
        case 1: return m_reader.m_spin_conserving_1e;
        case 2: return m_reader.m_spin_conserving_2e;
    }
    return true;
}

Hdf5IntegralReader::Hdf5IntegralReader(const FcidumpInfo& info, Hdf5IntegralReader::KeyNames names) :
        m_reader(info.m_fname), m_names(std::move(names)),
        m_indices_2e(m_reader, m_names.m_2e_inds),
        m_values_2e(m_reader.read_data<v_t<ham_comp_t>>(m_names.m_2e_values)),
        m_indices_1e(m_reader, m_names.m_1e_inds),
        m_values_1e(m_reader.read_data<v_t<ham_comp_t>>(m_names.m_1e_values)){
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

bool Hdf5IntegralReader::spin_conserving(uint_t /*iex*/) const {
    return true;
}

MolcasHdf5IntegralReader::MolcasHdf5IntegralReader(const FcidumpInfo& info) : Hdf5IntegralReader(info,
     {"CORE_ENERGY", "FOCK_INDEX", "FOCK_VALUES", "TWO_EL_INT_INDEX", "TWO_EL_INT_VALUES"}){}
