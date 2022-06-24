//
// Created by Robert J. Anderson on 03/06/2021.
//

#include "Shift.h"

Reweighter::Reweighter(const conf::Shift &opts, const NdFormat<c_ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_const_shift(wf_fmt.m_shape, opts.m_init),
        m_active("accumulating reweighting statistics", wf_fmt.m_nelement, "WF part"),
        m_histories(wf_fmt.m_nelement, std::queue<ham_comp_t>()),
        m_total(wf_fmt.m_shape, 1.0) {
}

void Reweighter::update(uint_t icycle, uint_t ipart, bool begin_cond, ham_comp_t av_shift) {
    if (m_opts.m_reweight.m_ncycle == 0) return;
    if (m_active[ipart].update(icycle, begin_cond)) {
        m_const_shift[ipart] = av_shift;
        log::info("setting constant shift for reweighting on WF part {} to current average shift {}",
                  ipart, m_const_shift[ipart]);
    }
}

void Reweighter::add(uint_t ipart, ham_comp_t shift, double tau) {
    if (m_opts.m_reweight.m_ncycle == 0) return;
    if (!m_active[ipart]) return;
    add_to_history(ipart, std::exp(tau * (m_const_shift[ipart] - shift)));
    ASSERT(product_chk());
    ASSERT(m_histories[ipart].size() <= m_opts.m_reweight.m_ncycle);
}

void Reweighter::add_to_history(uint_t ipart, const ham_comp_t &v) {
    auto &history = m_histories[ipart];
    history.push(v);
    m_total[ipart] *= v;
    if (history.size() == m_opts.m_reweight.m_ncycle) {
        auto oldest_factor = history.front();
        m_total[ipart] /= oldest_factor;
        history.pop();
    }
    // else we're still filling.
}

bool Reweighter::product_chk() const {
    for (uint_t ipart=0ul; ipart<m_const_shift.nelement(); ++ipart) {
        ham_comp_t total = 1.0;
        auto histories_copy = m_histories[ipart];
        while (!histories_copy.empty()) {
            total *= histories_copy.front();
            histories_copy.pop();
        }
        if (!fptol::numeric_equal(total, m_total[ipart])) return false;
    }
    return true;
}

Shift::Shift(const conf::Document &opts, const NdFormat<c_ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_nwalker_last_period(wf_fmt.m_shape, std::numeric_limits<wf_comp_t>::max()),
        m_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_avg_value_histories(wf_fmt.m_nelement, std::queue<ham_comp_t>()),
        m_avg_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_variable_mode("variable shift mode", wf_fmt.m_nelement, "WF part"),
        m_nwalker_target("nwalker_target", opts.m_propagator.m_nw_target),
        m_reweighter(opts.m_shift, wf_fmt){
    m_nwalker_last_period.zero();
    DEBUG_ASSERT_FALSE(m_variable_mode, "Shift should not initially be in variable mode");
}

const ham_comp_t &Shift::operator[](const uint_t &ipart) {
    return m_values[ipart];
}

void Shift::update(const Wavefunction &wf, const uint_t &icycle, const double &tau) {
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
            log::info("Variable shift triggered for WF part {}. Cycle {} nw: {}, cycle {} nw: {}",
                      ipart, icycle-1, wf.m_nwalker.m_reduced[ipart], icycle, nw);
            a = icycle % m_opts.m_shift.m_period;
        }
        if (is_period_cycle) a = m_opts.m_shift.m_period;

        if (variable_mode && a) {
            auto rate =  nw / m_nwalker_last_period[ipart];
            m_values[ipart] -= m_opts.m_shift.m_damp * std::log(std::abs(rate)) / (tau * a);
            if (m_opts.m_shift.m_reweight.enabled()) {
                bool reweight_begin_cond = icycle >= variable_mode.icycle_start() + m_opts.m_shift.m_reweight.m_delay;
                m_reweighter.update(icycle, ipart, reweight_begin_cond, get_average(ipart));
                m_reweighter.add(ipart, m_values[ipart], tau);
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
