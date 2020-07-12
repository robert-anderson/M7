//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"
#include "FciqmcCalculation.h"

Propagator::Propagator(FciqmcCalculation *fciqmc) :
    m_fciqmc(fciqmc), m_input(fciqmc->m_input),
    m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator),
    m_magnitude_logger(m_input),
    m_dst_det(m_ham->nsite()),
    m_aconn(m_dst_det), m_occ(m_dst_det), m_vac(m_dst_det),
    m_variable_shift(fciqmc->m_vary_shift),
    m_semi_stochastic(fciqmc->m_semi_stochastic)
    {
    m_tau = m_input.tau_initial;
    m_shift = m_input.shift_initial;
}

