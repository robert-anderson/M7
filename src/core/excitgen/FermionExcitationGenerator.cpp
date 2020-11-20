//
// Created by RJA on 20/11/2020.
//

#include "FermionExcitationGenerator.h"

FermionExcitationGenerator::FermionExcitationGenerator(const Hamiltonian *h, PRNG &prng, size_t nexcit) :
        ExcitationGenerator(h, prng), m_nexcit(nexcit),
        m_spin_conserving(nexcit==1 ? h->spin_conserving_1e() : h->spin_conserving_2e()){}

bool FermionExcitationGenerator::draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                                            const OccupiedOrbitals &occ, const VacantOrbitals &vac, defs::prob_t &prob,
                                            defs::ham_t &helem, conn::AsFermionOnv &anticonn) {
    return false;
}

bool FermionExcitationGenerator::draw(const views::FermiBosOnv &src_fonv, views::FermiBosOnv &dst_fonv,
                                            const OccupiedOrbitals &occ, const VacantOrbitals &vac, defs::prob_t &prob,
                                            defs::ham_t &helem, conn::AsFermiBosOnv &anticonn) {
    return draw(src_fonv.m_fonv, dst_fonv.m_fonv, occ, vac, prob, helem, anticonn.m_aconn);
}
