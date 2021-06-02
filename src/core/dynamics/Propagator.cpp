//
// Created by Robert John Anderson on 2020-02-11.
//


#include "Propagator.h"

Shift::Shift(const Options &opts, const NdFormat<defs::ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_nwalker_last_update(wf_fmt.shape(), std::numeric_limits<defs::wf_t>::max()),
        m_values(wf_fmt.shape(), opts.shift_initial),
        m_avg_values(wf_fmt.shape(), 0.0),
        m_const_shift(wf_fmt.shape(), opts.shift_initial),
        m_variable_mode("variable shift mode", wf_fmt.nelement(), "WF part"),
        m_reweighting_active("accumulating reweighting statistics", wf_fmt.nelement(), "WF part"),
        m_reweighting_factors(wf_fmt.nelement(), std::queue<defs::ham_t>()),
        m_total_reweighting(wf_fmt.shape(), 1.0),
        m_nwalker_target("nwalker_target", opts.nwalker_target){}

const defs::ham_comp_t &Shift::operator[](const size_t &ipart) {
    return m_values[ipart];
}


void Shift::evaluate_reweighting(const size_t& ipart, const size_t
&icycle, const double& tau) {
    if (m_opts.ncycle_reweight_lookback == 0) return;
    if (!m_reweighting_active[ipart]) return;
    ASSERT(m_variable_mode[ipart] && m_reweighting_active[ipart])

    auto this_factor = std::exp(tau * (m_const_shift[ipart] - m_values[ipart]));
    m_reweighting_factors[ipart].push(this_factor);
    m_total_reweighting[ipart] *= this_factor;

    if(m_reweighting_factors.size() == m_opts.ncycle_reweight_lookback){
        auto oldest_factor = m_reweighting_factors[ipart].front();
        m_total_reweighting[ipart] /= oldest_factor;
        m_reweighting_factors[ipart].pop();
    }
    // else we're still filling.
    ASSERT(consts::floats_nearly_equal(m_total_reweighting[ipart], product_reweighting_queue(ipart)))
    ASSERT(m_reweighting_factors.size() <= m_opts.ncycle_reweight_lookback)
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
            auto& reweighting_active = m_reweighting_active[ipart];
            bool begin_reweighting = m_opts.ncycle_reweight_lookback > 0;
            begin_reweighting &= (icycle >= variable_mode.icycle_start() + m_opts.ncycle_wait_reweight);
            // compute average shift from variable shift epoch up to now
            m_avg_values[ipart] += m_values[ipart];
            if(reweighting_active.update(icycle, begin_reweighting)){
                //m_const_shift[ipart] = m_avg_values[ipart] / (m_opts.ncycle_wait_reweight + 1);
                m_const_shift[ipart] = -2.9925653194551863;
                log::info("setting constant shift for reweighting on WF part {} to current average variable "
                          "shift from last {} cycles: C[{}] = {}", ipart, m_opts.ncycle_wait_reweight + 1,  ipart,
                          m_const_shift[ipart]);
            }
            evaluate_reweighting(ipart, icycle, tau);
        }
    }
    m_nwalker_last_update = wf.m_nwalker.m_reduced;
}


defs::ham_comp_t Shift::product_reweighting_queue(const size_t ipart) {
    defs::ham_comp_t total = 1.0;
    auto factors_copy = m_reweighting_factors[ipart];
    while(!factors_copy.empty()){
        total *= factors_copy.front();
        factors_copy.pop();
    }
    return total;
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