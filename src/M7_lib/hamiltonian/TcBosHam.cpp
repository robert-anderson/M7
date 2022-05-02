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

defs::ham_t TcBosHam::get_coeff_0033(size_t a, size_t b, size_t c,
                                     size_t i, size_t j, size_t k) const {
    const int ia=i, ib=b, ic=c, ii=i, ij=j, ik=k; // convert to int
    return three_body_coeff(&ia, &ib, &ic, &ii, &ij, &ik);
}

defs::ham_t TcBosHam::get_element_0000(const field::BosOnv &onv) const {
    auto element = GeneralBosHam::get_element_0000(onv);
    // @todo stub requires sum over triples
}

defs::ham_t TcBosHam::get_element_0011(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    auto element = GeneralBosHam::get_element_0011(onv, conn);
    // TODO stub
}

defs::ham_t TcBosHam::get_element_0022(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    auto element = GeneralBosHam::get_element_0022(onv, conn);
    for(size_t imode=0; imode < onv.m_nmode; ++imode) {
        // no need to check operators like in Fermi case
        element += get_coeff_0033(conn.m_cre.get_imode(0), conn.m_cre.get_imode(1), imode,
                                  conn.m_ann.get_imode(0), conn.m_ann.get_imode(1), imode);
        // CHECK not sure if need to add other repeated versions, like if two m_cre objects are repeated
    }
}

defs::ham_t TcBosHam::get_element_0033(const field::BosOnv &onv, const conn::BosOnv &conn) const {
    auto element = get_coeff_0033(conn.m_cre.get_imode(0), conn.m_cre.get_imode(1), conn.m_cre.get_imode(2),
                                  conn.m_ann.get_imode(0), conn.m_ann.get_imode(1), conn.m_ann.get_imode(2));
    element *= conn.occ_fac(onv);
}

