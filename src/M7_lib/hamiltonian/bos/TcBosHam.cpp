/**
 * @file TcBosHam.cpp
 * @author jph
 * @brief transcorrelated Boson Hamiltonian implementation file
 * @date 2022-04-22
 *
 */

#include "TcBosHam.h"

// note that in the Boson case we need not (cannot) sum over bit triples as we
// can have a higher occupation than one for a given mode

defs::ham_t TcBosHam::get_coeff_0033(uint_t a, uint_t b, uint_t c, uint_t i,
                                     uint_t j, uint_t k) const {
    // symmetrise
    defs::ham_t coeff = get_lmat_coeff(a, b, c, i, j, k)
                      + get_lmat_coeff(a, b, c, j, k, i)
                      + get_lmat_coeff(a, b, c, k, i, j)
                      + get_lmat_coeff(a, b, c, j, i, k)
                      + get_lmat_coeff(a, b, c, i, k, j)
                      + get_lmat_coeff(a, b, c, k, j, i);
    return coeff;
}

defs::ham_t TcBosHam::get_element_0000(const field::BosOnv &onv) const {
    auto element = GeneralBosHam::get_element_0000(onv);
    for (uint_t imode = 0; imode < onv.m_basis.m_nmode; ++imode) {
        for (uint_t jmode = 0; jmode < onv.m_basis.m_nmode; ++jmode) {
            for (uint_t kmode = 0; kmode < onv.m_basis.m_nmode; ++kmode) {
                element +=
                    get_coeff_0033(imode, jmode, kmode, imode, jmode, kmode);
            }
        }
    }
    return element;
}

defs::ham_t TcBosHam::get_element_0011(const field::BosOnv &onv,
                                       const conn::BosOnv &conn) const {
    auto element = GeneralBosHam::get_element_0011(onv, conn);
    for (uint_t imode = 0; imode < onv.m_basis.m_nmode; ++imode) {
        for (uint_t jmode = 0; jmode < onv.m_basis.m_nmode; ++jmode) {
            // CHECK am I double counting here?
            BosOps com(2);
            com.set(imode, jmode);
            element += onv.occ_fac(conn, com) *
                       get_coeff_0033(conn.m_cre.get_imode(0), imode, jmode,
                                      conn.m_ann.get_imode(0), imode, jmode);
            // CHECK are other permutations necessary to add?
        }
    }
    return element;
}

defs::ham_t TcBosHam::get_element_0022(const field::BosOnv &onv,
                                       const conn::BosOnv &conn) const {
    auto element = GeneralBosHam::get_element_0022(onv, conn);
    for (uint_t imode = 0; imode < onv.m_basis.m_nmode; ++imode) {
        // no need to check operators like in Fermi case
        BosOps com(1);
        com.set(imode);
        element += onv.occ_fac(conn, com) *
                   get_coeff_0033(
                       conn.m_cre.get_imode(0), conn.m_cre.get_imode(1), imode,
                       conn.m_ann.get_imode(0), conn.m_ann.get_imode(1), imode);
        // CHECK not sure if need to add other repeated versions, like if two
        // m_cre objects are repeated
    }
    return element;
}

defs::ham_t TcBosHam::get_element_0033(const field::BosOnv &onv,
                                       const conn::BosOnv &conn) const {
    auto element =
        get_coeff_0033(conn.m_cre.get_imode(0), conn.m_cre.get_imode(1),
                       conn.m_cre.get_imode(2), conn.m_ann.get_imode(0),
                       conn.m_ann.get_imode(1), conn.m_ann.get_imode(2));
    element *= onv.occ_fac(conn);
    return element;
}
