//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSHAM_H
#define M7_BOSHAM_H

#include <src/core/integrals/BosonCoeffs_1.h>
#include <src/core/integrals/BosonCoeffs_2.h>
#include <src/core/io/BosdumpFileReader.h>
#include "src/core/connection/Connections.h"
#include "src/core/parallel/SharedArray.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct BosHam {
    const size_t m_nmode, m_nboson;
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;
    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    BosHam(const BosdumpHeader& header);

    BosHam(const std::string& fname): BosHam(BosdumpHeader(fname)){}

    defs::ham_t get_element(const field::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const;

    defs::ham_t get_element(const field::BosOnv &src, const conn::BosOnv& conn) const;

    size_t nci() const;

    void log_data() const;

};


#endif //M7_BOSHAM_H
