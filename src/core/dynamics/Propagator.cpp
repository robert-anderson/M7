//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"
#include "FciqmcCalculation.h"

Propagator::Propagator(FciqmcCalculation *fciqmc) :
    m_fciqmc(fciqmc), m_input(fciqmc->m_input),
    m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator) {
    m_tau = m_input.tau_initial;
    m_shift = m_input.shift_initial;

}

