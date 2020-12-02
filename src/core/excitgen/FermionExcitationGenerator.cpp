//
// Created by RJA on 20/11/2020.
//

#include "FermionExcitationGenerator.h"

FermionExcitationGenerator::FermionExcitationGenerator(const Hamiltonian<> *h, PRNG &prng, size_t nexcit) :
        ExcitationGenerator(h, prng), m_nexcit(nexcit),
        m_spin_conserving(nexcit==1 ? h->spin_conserving_1e() : h->spin_conserving_2e()){}

bool FermionExcitationGenerator::draw(const views::Onv<0> &src_onv, views::Onv<0> &dst_onv,
                                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                                      defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<0> &anticonn) {
    return false;
}

bool FermionExcitationGenerator::draw(const views::Onv<1> &src_onv, views::Onv<1> &dst_onv,
                                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                                      defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<1> &anticonn) {
    return draw(src_onv.m_fonv, dst_onv.m_fonv, occs, vacs, prob, helem, anticonn);
}
