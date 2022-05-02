//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALBOSHAM_H
#define M7_GENERALBOSHAM_H

#include <M7_lib/config/Hamiltonian.h>

#include "BosHam.h"

struct GeneralBosHam : BosHam {
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;

    GeneralBosHam(const BosdumpHeader &header);

    GeneralBosHam(const fciqmc_config::BosonHamiltonian &opts) : GeneralBosHam(BosdumpHeader(opts.m_bosdump.m_path)){}

    defs::ham_t get_coeff_0011(size_t i, size_t j) const override;

    defs::ham_t get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const override;

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override;

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

};

#endif //M7_GENERALBOSHAM_H

// todo remove comments below
// const, const ref, method const
// move semantics e.g. && -> might solve copy problem
// smart pointers
// garbage collection
// smart pointers
// **std::unique_ptr** (for polymorphic obj, singular ownership)
// std::share_ptr (not useful here? but useful in genreal, shared ownership)
// use standard template library e.g. std::vector
// virtual methods
