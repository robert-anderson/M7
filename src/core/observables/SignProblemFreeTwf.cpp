//
// Created by jhalson on 28/04/2021.
//

#include "SignProblemFreeTwf.h"

SpfTwfBase::SpfTwfBase(const Hamiltonian &ham, size_t npart, size_t nsite) :
        m_ham(ham),  m_excit_iters(ham), m_numerator(npart, 0.0), m_denominator(npart, 0.0),
        m_numerator_total(npart, 0.0), m_denominator_total(npart, 0.0) {
}

void SpfTwfBase::reduce() {
    mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
    mpi::all_sum(m_denominator.data(), m_denominator_total.data(), m_denominator.size());
    m_numerator.assign(m_numerator.size(), 0.0);
    m_denominator.assign(m_denominator.size(), 0.0);
}
