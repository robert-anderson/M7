//
// Created by RJA on 20/11/2020.
//

#include "FermionExcitationGenerator.h"

FermionExcitationGenerator::FermionExcitationGenerator(const Hamiltonian *h, PRNG &prng, size_t nexcit) :
        ExcitationGenerator(h, prng), m_nexcit(nexcit),
        m_spin_conserving(nexcit==1 ? h->m_frm.spin_conserving_1e() : h->m_frm.spin_conserving_2e()){}