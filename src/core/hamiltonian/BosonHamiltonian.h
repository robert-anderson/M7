//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONHAMILTONIAN_H
#define M7_BOSONHAMILTONIAN_H


#include <src/core/field/Fields.h>

struct BosonHamiltonian {
    const size_t m_nboson_cutoff, m_nmode;
    const defs::ham_t m_omega;

    BosonHamiltonian(size_t nmode, size_t nboson_cutoff, defs::ham_t omega);

    defs::ham_t get_element(const fields::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const fields::BosOnv &onv) const;
};


#endif //M7_BOSONHAMILTONIAN_H
