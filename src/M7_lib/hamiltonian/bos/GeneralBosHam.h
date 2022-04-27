//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALBOSHAM_H
#define M7_GENERALBOSHAM_H

#include "M7_lib/config/Hamiltonian.h"

#include "M7_lib/hamiltonian/bos/BosHam.h"

struct GeneralBosHam : BosHam {
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;

    GeneralBosHam(const BosdumpHeader &header, size_t occ_cutoff);

    GeneralBosHam(const conf::BosHam &opts);

    defs::ham_t get_coeff_0011(size_t i, size_t j) const override;

    defs::ham_t get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const override;

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override;

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

};

#endif //M7_GENERALBOSHAM_H
