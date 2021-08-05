//
// Created by rja on 18/05/2020.
//

#include "MagnitudeLogger.h"

MagnitudeLogger::MagnitudeLogger(defs::ham_comp_t max_bloom, size_t ndraw_min, size_t nexlvl, bool static_tau,
                                 bool static_probs, double tau_min, double tau_max, double prob_min, size_t period) :
        m_max_bloom(max_bloom), m_ndraw_min(ndraw_min), m_static_tau(static_tau), m_static_probs(static_probs),
        m_tau_min(tau_min), m_tau_max(tau_max), m_prob_min(prob_min), m_period(period), m_ndraw({nexlvl}),
        m_gamma({nexlvl}), m_new_probs(nexlvl){
    log::info("Initializing magnitude logger with max_bloom {} for {} excitation levels", max_bloom, nexlvl);
    log::info("Dynamic tau optimization: {}", !m_static_tau);
    if (!m_static_tau)
        log::info("Dynamic tau to be kept above {} and below {}", m_tau_min, m_tau_max);
    log::info("Dynamic excitation level probability optimization: {}", !m_static_probs);
    if (!m_static_probs)
        log::info("Dynamic excitation level probabilities to be kept above {}", m_prob_min);
}

void MagnitudeLogger::log(const size_t &iexlvl, const defs::ham_comp_t &helem, const defs::prob_t &prob) {
    DEBUG_ASSERT_NE(prob, 0.0, "null draw should never be logged");
    auto& hi = m_gamma.m_local[iexlvl];
    auto mag = helem / prob;
    if (mag > hi) hi = mag;
    ++m_ndraw.m_local[iexlvl];
}

void MagnitudeLogger::update_tau(double &tau, const defs::ham_comp_t &gamma_sum) {
    if (m_static_tau) return;
    tau = m_max_bloom / gamma_sum;
    if (tau < m_tau_min) tau = m_tau_min;
    if (tau > m_tau_max) tau = m_tau_max;
}

void MagnitudeLogger::update(size_t icycle, double &tau) {
    if (!icycle || icycle%m_period) return;
    m_gamma.all_max();
    update_tau(tau, m_gamma.m_reduced.sum());
}

void MagnitudeLogger::update(size_t icycle, double &tau, ExcitGenGroup &excit_gens) {
    if (!icycle || icycle%m_period) return;
    m_gamma.all_max();
    auto gamma_sum = m_gamma.m_reduced.sum();
    update_tau(tau, gamma_sum);
    if (!m_static_probs){
        m_new_probs.assign(m_gamma.m_reduced.dbegin(), m_gamma.m_reduced.dend());
        for (auto& new_prob : m_new_probs) new_prob /= gamma_sum;
        prob_utils::rectify(m_new_probs, m_prob_min);
        excit_gens.set_probs(m_new_probs);
        m_gamma.m_local.zero();
    }
}
