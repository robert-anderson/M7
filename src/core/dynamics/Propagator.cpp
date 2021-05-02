//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"

Shift::Shift(const Options &opts, const NdFormat<defs::ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_nwalker_last_update(wf_fmt.shape(), std::numeric_limits<defs::wf_t>::max()),
        m_values(wf_fmt.shape(), opts.shift_initial),
        m_variable_mode("variable shift mode", wf_fmt.nelement(), "WF part"),
        m_nwalker_target("nwalker_target", opts.nwalker_target){}

const defs::ham_comp_t &Shift::operator[](const size_t &ipart) {
    return m_values[ipart];
}

void Shift::update(const Wavefunction &wf, const size_t &icycle, const double &tau) {
    if (m_nwalker_target.read()) m_variable_mode.terminate(icycle);

    if (icycle % m_opts.shift_update_period) return;

    if (m_nwalker_last_update[0]==std::numeric_limits<defs::wf_t>::max()){
        m_nwalker_last_update = wf.m_nwalker.m_reduced;
        return;
    }

    for (size_t ipart=0ul; ipart < m_variable_mode.nelement(); ++ipart){
        auto& variable_mode = m_variable_mode[ipart];
        variable_mode.update(icycle, wf.m_nwalker.m_reduced[ipart] >= m_nwalker_target);
        if (variable_mode) {
            auto rate = wf.m_nwalker.m_reduced[ipart] / m_nwalker_last_update[ipart];
            m_values[ipart] -= m_opts.shift_damp * consts::real_log(rate) / (tau * m_opts.shift_update_period);
        }
    }
    m_nwalker_last_update = wf.m_nwalker.m_reduced;
}



//Propagator::Propagator() :
//        m_ham(fciqmc->m_ham), m_rank_allocator(fciqmc->m_rank_allocator),
//        m_magnitude_logger(m_input, m_ham->nsite(), m_ham->nelec()),
//        m_dst_det(m_ham->nsite()),
//        m_aconn(m_dst_det), m_occ(m_dst_det), m_vac(m_dst_det),
//        m_variable_shift(fciqmc->m_vary_shift),
//        m_semi_stochastic(fciqmc->m_semi_stochastic)
//{
//    m_values = m_input.shift_initial;
//    std::cout << "Propagator base initialized with shift (relative to reference determinant energy): " << m_values << std::endl;
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
//    m_values = m_input.shift_initial;
//    std::cout << "Propagator base initialized with shift (relative to reference determinant energy): " << m_values << std::endl;
//}

void Propagator::update(const size_t& icycle, const Wavefunction& wf) {
    //m_magnitude_logger.synchronize(icycle);
    m_shift.update(wf, icycle, tau());
}