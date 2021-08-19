//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONHAMILTONIAN_H
#define M7_BOSONHAMILTONIAN_H

#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"

struct BosonHamiltonian {
    const size_t m_nboson_max, m_nmode;
    const defs::ham_t m_omega;

    BosonHamiltonian(size_t nmode, size_t nboson_max, std::string fname);

    defs::ham_t get_element(const field::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const;

    defs::ham_t get_element(const field::BosOnv &onv, const conn::BosOnv& conn) const {
        if (conn.size()) return 0.0;
        return get_element(onv);
    }

    size_t nci() const;
};


#endif //M7_BOSONHAMILTONIAN_H
