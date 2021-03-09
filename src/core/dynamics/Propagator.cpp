//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"

//Propagator::Propagator() :
//        m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator),
//        m_magnitude_logger(m_input, m_ham->nsite(), m_ham->nelec()),
//        m_dst_det(m_ham->nsite()),
//        m_aconn(m_dst_det), m_occ(m_dst_det), m_vac(m_dst_det),
//        m_variable_shift(fciqmc->m_vary_shift),
//        m_semi_stochastic(fciqmc->m_semi_stochastic)
//{
//    m_shift = m_input.shift_initial;
//    std::cout << "Propagator base initialized with shift (relative to reference determinant energy): " << m_shift << std::endl;
//}

#include "Propagator.h"
#include "FciqmcCalculation.h"

//Propagator::Propagator(FciqmcCalculation *fciqmc) :
//    m_fciqmc(fciqmc), m_input(fciqmc->m_input),
//    m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator),
//    m_magnitude_logger(m_input, m_ham->nsite(), m_ham->nelec()),
//    m_dst_det(m_ham->nsite()),
//    m_aconn(m_dst_det), m_occ(m_dst_det), m_vac(m_dst_det),
//    m_variable_shift(fciqmc->m_vary_shift),
//    m_semi_stochastic(fciqmc->m_semi_stochastic)
//    {
//    m_shift = m_input.shift_initial;
//    std::cout << "Propagator base initialized with shift (relative to reference determinant energy): " << m_shift << std::endl;
//}

void Propagator::update(const size_t& icycle, const Wavefunction& wf) {
    //m_magnitude_logger.synchronize(icycle);
    m_variable_shift.update(icycle, wf.m_nwalker.reduced(0, 0) >= m_opts.nwalker_target);
    if (icycle % m_opts.shift_update_period) return;
//    if (m_variable_shift.update(icycle, wf.m_nwalker.reduced() >= m_opts.nwalker_target)) {
//        /*
//         * "jump" the shift to the projected energy estimation at the onset of
//         * the variable shift epoch
//         */
//        m_shift = wf.refref_proj_energy();
//    }
    else if (m_variable_shift) {
        auto rate = 1.0+wf.m_delta_nwalker.reduced(0, 0)/wf.m_nwalker.reduced(0, 0);
        m_shift -= m_opts.shift_damp * consts::real_log(rate) / (tau()*m_opts.shift_update_period);
    }
}
