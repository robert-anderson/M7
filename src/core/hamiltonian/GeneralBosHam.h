//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALBOSHAM_H
#define M7_GENERALBOSHAM_H

#include "BosHam.h"

#if 0
class GeneralBosHam : BosHam {
    const size_t m_nmode, m_nboson;
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;
    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    BosHam(const BosdumpHeader &header);

    BosHam(const std::string &fname) : BosHam(BosdumpHeader(fname)) {}

    defs::ham_t get_element(const field::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const;

    defs::ham_t get_element(const field::BosOnv &src, const conn::BosOnv &conn) const;

    size_t nci() const;

    void log_data() const;

};

#endif

#endif //M7_GENERALBOSHAM_H
