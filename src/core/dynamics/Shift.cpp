//
// Created by rja on 03/06/2021.
//

#include "Shift.h"

Reweighter::Reweighter(const fciqmc_config::Shift &opts, const NdFormat<defs::ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_const_shift(wf_fmt.m_shape, opts.m_init),
        m_active("accumulating reweighting statistics", wf_fmt.m_nelement, "WF part"),
        m_histories(wf_fmt.m_nelement, std::queue<defs::ham_comp_t>()),
        m_total(wf_fmt.m_shape, 1.0) {
}

void Reweighter::update(size_t icycle, size_t ipart, bool begin_cond, defs::ham_comp_t av_shift) {
    if (m_opts.m_reweight.m_ncycle == 0) return;
    if (m_active[ipart].update(icycle, begin_cond)) {
        m_const_shift[ipart] = av_shift;
        log::info("setting constant shift for reweighting on WF part {} to current average shift {}",
                  ipart, m_const_shift[ipart]);
    }
}

void Reweighter::add(size_t ipart, defs::ham_comp_t shift, size_t icycle, double tau) {
    if (m_opts.m_reweight.m_ncycle == 0) return;
    if (!m_active[ipart]) return;
    add_to_history(ipart, std::exp(tau * (m_const_shift[ipart] - shift)));
    ASSERT(product_chk());
    ASSERT(m_histories[ipart].size() <= m_opts.m_reweight.m_ncycle);
}

void Reweighter::add_to_history(size_t ipart, const defs::ham_comp_t &v) {
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
    for (size_t ipart=0ul; ipart<m_const_shift.nelement(); ++ipart) {
        defs::ham_comp_t total = 1.0;
        auto histories_copy = m_histories[ipart];
        while (!histories_copy.empty()) {
            total *= histories_copy.front();
            histories_copy.pop();
        }
        if (!consts::floats_nearly_equal(total/m_total[ipart], 1.0, 1e-12))
            return false;
    }
    return true;
}

Shift::Shift(const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf> &wf_fmt) :
        m_opts(opts),
        m_nwalker_last_update(wf_fmt.m_shape, std::numeric_limits<defs::wf_t>::max()),
        m_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_avg_value_histories(wf_fmt.m_nelement, std::queue<defs::ham_comp_t>()),
        m_avg_values(wf_fmt.m_shape, opts.m_shift.m_init),
        m_variable_mode("variable shift mode", wf_fmt.m_nelement, "WF part"),
        m_nwalker_target("nwalker_target", opts.m_propagator.m_nw_target),
        m_reweighter(opts.m_shift, wf_fmt){}

const defs::ham_comp_t &Shift::operator[](const size_t &ipart) {
    return m_values[ipart];
}

void Shift::update(const Wavefunction &wf, const size_t &icycle, const double &tau) {
    if (m_nwalker_target.read()) m_variable_mode.terminate(icycle);

    if (icycle % m_opts.m_shift.m_period) return;

    /**
     * we are computing S(icycle), but have the walker number at icycle-1 Nw(icycle-1)
     * to compute S(icycle), we need Nw(icycle)
     *
     * this is constructed by Nw(icycle-1)+deltaNw(icycle-1)
     *
     * the ratio between Nw(icycle) and Nw(icycle - period) is taken to determine the new shift value
     *
     * Nw(icycle) is stored in m_nwalker_last_update
     *
     */
     /**
      * on the 0-th cycle this makes no difference, since Nw(-period) is taken to be zero
      */
    if (m_nwalker_last_update[0]==std::numeric_limits<defs::wf_t>::max()){
        m_nwalker_last_update = wf.m_nwalker.m_reduced;
        return;
    }

    for (size_t ipart=0ul; ipart < m_variable_mode.nelement(); ++ipart){
        auto& variable_mode = m_variable_mode[ipart];
        std::cout << wf.m_nwalker.m_reduced[ipart] << " " << m_nwalker_target << std::endl;
        variable_mode.update(icycle, wf.m_nwalker.m_reduced[ipart] >= m_nwalker_target);
        if (variable_mode) {
            auto next_nw = wf.m_nwalker.m_reduced[ipart] + wf.m_delta_nwalker.m_reduced[ipart];
            auto rate = next_nw / m_nwalker_last_update[ipart];
            m_values[ipart] -= m_opts.m_shift.m_damp * consts::real_log(rate) / (tau * m_opts.m_shift.m_period);
            if (m_opts.m_shift.m_reweight) {
                bool reweight_begin_cond = icycle >= variable_mode.icycle_start() + m_opts.m_shift.m_reweight.m_delay;
                m_reweighter.update(icycle, ipart, reweight_begin_cond, get_average(ipart));
                m_reweighter.add(ipart, m_values[ipart], icycle, tau);
            }

        }
        add_to_average();
        m_nwalker_last_update[ipart] = wf.m_nwalker.m_reduced[ipart] + wf.m_delta_nwalker.m_reduced[ipart];
    }
}

void Shift::add_to_average() {
    for (size_t ipart=0ul; ipart<m_values.nelement(); ++ipart) {
        m_avg_value_histories[ipart].push(m_values[ipart]);
        if (m_avg_value_histories[ipart].size() > m_opts.m_shift.m_ncycle_av) {
            m_avg_values[ipart] -= m_avg_value_histories[ipart].front();
            m_avg_value_histories[ipart].pop();
        }
        ASSERT(m_avg_value_histories[ipart].size() <= m_opts.m_shift.m_ncycle_av)
        m_avg_values[ipart] += m_values[ipart];
    }
}

defs::ham_comp_t Shift::get_average(const size_t &ipart) const {
    return m_avg_values[ipart]/m_avg_value_histories[ipart].size();
}
