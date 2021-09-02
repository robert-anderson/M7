//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONHAMILTONIAN_H
#define M7_BOSONHAMILTONIAN_H

#include <src/core/integrals/BosonCoeffs.h>
#include "src/core/connection/Connections.h"
#include "src/core/parallel/SharedArray.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct BosonHamiltonian {
    const size_t m_nboson_max, m_nmode;
    BosonCoeffs m_coeffs;
    ham_data::TermContribs m_contribs_0011;

    BosonHamiltonian(size_t nmode, size_t nboson_max, std::string fname);

    defs::ham_t get_element(const field::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const;

    defs::ham_t get_element(const field::BosOnv &onv, const conn::BosOnv& conn) const;

    size_t nci() const;

    void log_data() const;
};


#endif //M7_BOSONHAMILTONIAN_H
