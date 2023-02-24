//
// Created by Robert J. Anderson on 03/06/2021.
//

#include "Shift.h"
#include "Wavefunction.h"
#include "Reference.h"

Shift::Shift(const conf::Document &opts, const NdFormat<c_ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_nwalker_last_period(wf_fmt.m_shape, std::numeric_limits<wf_comp_t>::max()),
        m_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_avg_value_histories(wf_fmt.m_nelement, std::queue<ham_comp_t>()),
        m_avg_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_variable_mode("variable shift mode", wf_fmt.m_nelement, "WF part"),
        m_nwalker_target("nwalker_target", opts.m_propagator.m_nw_target){
    if (m_opts.m_shift.m_target_damp)
        logging::info("Using targeted shift damping term (growth relative to target population)");
    m_nwalker_last_period.zero();
    DEBUG_ASSERT_FALSE(m_variable_mode, "Shift should not initially be in variable mode");
}

const ham_comp_t &Shift::operator[](uint_t ipart) {
    return m_values[ipart];
}

void Shift::update(const wf::Fci &wf, const wf::References& refs, uint_t icycle, double tau) {
    if (m_nwalker_target.read()) m_variable_mode.terminate(icycle);
    const bool is_period_cycle = !(icycle % m_opts.m_shift.m_period);

    for (uint_t ipart=0ul; ipart < m_variable_mode.nelement(); ++ipart){
        /*
         * at the beginning of cycle i - where this update is performed, Nw_i is not available directly since the loop
         * over the current occupied list has yet to be performed. Nw_i is required so compute S_i, so we must get it by
         * adding the difference in Nw due to the application of cycle i-1 propagator.
         */
        auto nw = wf.m_nwalker.m_reduced[ipart] + wf.m_delta_nwalker.m_reduced[ipart];
        auto& variable_mode = m_variable_mode[ipart];
        /*
         * number of cycles since last update
         */
        uint_t a = 0ul;
        if (variable_mode.update(icycle, nw >= m_nwalker_target)) {
            logging::info("Variable shift triggered for WF part {}. Cycle {} nw: {}, cycle {} nw: {}",
                      ipart, icycle-1, wf.m_nwalker.m_reduced[ipart], icycle, nw);
            a = icycle % m_opts.m_shift.m_period;
        }
        if (is_period_cycle) a = m_opts.m_shift.m_period;

        if (variable_mode && a) {
            if (m_opts.m_shift.m_fix_ref_weight.m_value) {
                m_values[ipart] = refs[ipart].proj_energy_num()/refs[ipart].weight();
            }
            else {
                auto rate = nw / m_nwalker_last_period[ipart];
                m_values[ipart] -= m_opts.m_shift.m_damp * std::log(std::abs(rate)) / (tau * a);
                if (m_opts.m_shift.m_target_damp) {
                    rate = nw / m_nwalker_target;
                    m_values[ipart] -= m_opts.m_shift.m_target_damp * std::log(std::abs(rate)) / (tau * a);
                }
            }
        }
        add_to_average();
    }

    if (is_period_cycle){
        m_nwalker_last_period = wf.m_nwalker.m_reduced;
        m_nwalker_last_period += wf.m_delta_nwalker.m_reduced;
    }
}

void Shift::add_to_average() {
    for (uint_t ipart=0ul; ipart<m_values.nelement(); ++ipart) {
        m_avg_value_histories[ipart].push(m_values[ipart]);
        if (m_avg_value_histories[ipart].size() > m_opts.m_shift.m_ncycle_av) {
            m_avg_values[ipart] -= m_avg_value_histories[ipart].front();
            m_avg_value_histories[ipart].pop();
        }
        ASSERT(m_avg_value_histories[ipart].size() <= m_opts.m_shift.m_ncycle_av)
        m_avg_values[ipart] += m_values[ipart];
    }
}

ham_comp_t Shift::get_average(const uint_t &ipart) const {
    return m_avg_values[ipart]/m_avg_value_histories[ipart].size();
}
