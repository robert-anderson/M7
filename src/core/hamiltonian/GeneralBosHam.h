//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALBOSHAM_H
#define M7_GENERALBOSHAM_H

#include <config/Hamiltonian.h>

#include "BosHam.h"

struct GeneralBosHam : BosHam {
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;

    GeneralBosHam(const BosdumpHeader &header);

    GeneralBosHam(const fciqmc_config::BosonHamiltonian &opts) : GeneralBosHam(BosdumpHeader(opts.m_bosdump.m_path)){}

    defs::ham_t get_coeff_0011(const size_t &i, const size_t &j) const override;

    defs::ham_t get_coeff_0022(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override;

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

};

#endif //M7_GENERALBOSHAM_H
