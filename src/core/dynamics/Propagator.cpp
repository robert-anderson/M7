//
// Created by Robert John Anderson on 2020-02-11.
//

#include "Propagator.h"
#include "FciqmcCalculation.h"

Propagator::Propagator(FciqmcCalculation *fciqmc) :
    m_fciqmc(fciqmc), m_input(fciqmc->m_input),
    m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator),
    m_magnitude_logger(m_input, m_ham->nsite(), m_ham->nelec()),
    m_dst_det(m_ham->nsite()),
    m_aconn(m_dst_det), m_occ(m_dst_det), m_vac(m_dst_det),
    m_variable_shift(fciqmc->m_vary_shift),
    m_semi_stochastic(fciqmc->m_semi_stochastic)
    {
    m_shift = m_input.shift_initial;
    std::cout << "Propagator initialized with shift (relative to reference determinant energy): " << m_shift << std::endl;
}

void Propagator::update(const size_t& icycle, defs::wf_comp_t nwalker, defs::wf_comp_t nwalker_growth) {
    m_magnitude_logger.synchronize(icycle);
    if (icycle % m_input.shift_update_period) return;
    if (m_variable_shift.update(icycle, nwalker >= m_input.nwalker_target)) {
        /*
         * "jump" the shift to the projected energy estimation at the onset of
         * the variable shift epoch
         */
        m_shift = m_fciqmc->m_wf.ref_proj_energy();
    }
    else if (m_variable_shift) {
        m_shift -= m_input.shift_damp * consts::real_log(nwalker_growth) / tau();
    }
}

