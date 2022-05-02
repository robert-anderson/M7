/**
 * @file TcBosHam.h
 * @author jph
 * @brief transcorrelated bosonic Hamiltonian struct
 * @date 2022-04-20
 *
 */

#ifndef M7_TCBOSHAM_H
#define M7_TCBOSHAM_H

#include "GeneralBosHam.h"

// @todo this will actually come from an external module (tchint)
extern float three_body_coeff(const int* ida, const int* idb, const int* idc, \
                              const int* idi, const int* idj, const int* idk);

struct TcBosHam : GeneralBosHam {
    // as in the Fermion case, need to write a get_coeff_0033 method and also
    // need to rewrite the get_element methods as they get the 0033 coeffs
    // "folded in"
    // note we do not ned to rewrite the get_coeff methods as we have a K+U
    // FCIDUMP file
    explicit TcBosHam(const fciqmc_config::BosonHamiltonian &opts) : GeneralBosHam(opts) {}

    // defs::ham_t get_coeff_0011(size_t i, size_t j) const override;

    // defs::ham_t get_coeff_0022(size_t i, size_t j, size_t k, size_t l) const override;

    /**
     * @brief Get the coeff of type 0033 (Bosonic "triple excitation")
     */
    defs::ham_t get_coeff_0033(size_t a, size_t b, size_t c,
                               size_t i, size_t j, size_t k) const;

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override;

    defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0033(const field::BosOnv &onv, const conn::BosOnv &conn) const;

};

#endif // M7_TCBOSHAM_H