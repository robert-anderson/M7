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
#include "M7_lib/hamiltonian/TcHam.h"

struct TcBosHam : TcHam, GeneralBosHam {

    TcBosHam(const BosdumpHeader &header, uint_t occ_cutoff) : TcHam(), GeneralBosHam(header, occ_cutoff) {};

    // as in the Fermion case, need to write a get_coeff_0033 method and also
    // need to rewrite the get_element methods as they get the 0033 coeffs
    // "folded in"
    // note we do not ned to rewrite the get_coeff methods as we have a K+U
    // FCIDUMP file
    explicit TcBosHam(opt_pair_t opts) : TcHam(), GeneralBosHam(opts) {}

    /**
     * @brief Get the coeff of type 0033 (Bosonic "triple excitation")
     */
    defs::ham_t get_coeff_0033(uint_t a, uint_t b, uint_t c, uint_t i, uint_t j,
                               uint_t k) const;

    defs::ham_t get_element_0000(const field::BosOnv &onv) const override;

    defs::ham_t get_element_0011(const field::BosOnv &onv,
                                 const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0022(const field::BosOnv &onv,
                                 const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0033(const field::BosOnv &onv,
                                 const conn::BosOnv &conn) const;
};

#endif  // M7_TCBOSHAM_H