//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"
#include "FciqmcCalculation.h"

Propagator::Propagator(FciqmcCalculation *fciqmc):
    m_fciqmc(fciqmc), m_input(fciqmc->m_input), 
    m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator) {
    m_tau = m_input.tau_initial;
    m_shift = m_input.shift_initial;
}

void Propagator::diagonal(
    const NumericElement<defs::ham_comp_t> &hdiag, NumericElement<defs::ham_t> &weight, 
    defs::ham_comp_t &delta_square_norm, defs::ham_comp_t &delta_nw) const {
    auto tmp = *weight;
    weight *= (1.0 - (*hdiag - m_shift) * m_tau);
    delta_square_norm += std::pow(std::abs(*weight), 2) - std::pow(std::abs(tmp), 2);
    delta_nw += std::abs(*weight) - std::abs(tmp);
}
