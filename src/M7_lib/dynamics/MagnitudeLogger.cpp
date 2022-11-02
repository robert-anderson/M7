//
// Created by Robert J. Anderson on 18/05/2020.
//

#include "M7_lib/util/Probability.h"
#include "MagnitudeLogger.h"

MagnitudeLogger::MagnitudeLogger(ham_comp_t max_bloom, uint_t ndraw_min, uint_t nexcase, bool static_tau,
                                 bool static_probs, double tau_min, double tau_max, double prob_min, uint_t period) :
        m_max_bloom(max_bloom), m_ndraw_min(ndraw_min), m_static_tau(static_tau), m_static_probs(static_probs),
        m_tau_min(tau_min), m_tau_max(tau_max), m_prob_min(prob_min), m_period(period), m_ndraw({nexcase}),
        m_gamma({nexcase}), m_new_probs(nexcase){
    logging::info("Initializing magnitude logger with max_bloom {} for {} excitation levels", max_bloom, nexcase);
    logging::info("Dynamic tau optimization: {}", !m_static_tau);
    if (!m_static_tau)
        logging::info("Dynamic tau to be kept above {} and below {}", m_tau_min, m_tau_max);
    logging::info("Dynamic excitation level probability optimization: {}", !m_static_probs);
    if (!m_static_probs)
        logging::info("Dynamic excitation level probabilities to be kept above {}", m_prob_min);
}

void MagnitudeLogger::log(uint_t icase, const ham_t &helem, const prob_t &prob) {
    DEBUG_ASSERT_NE(prob, 0.0, "null draw should never be logged");
    auto& hi = m_gamma.m_local[icase];
    auto mag = std::abs(helem) / prob;
    if (mag > hi) hi = mag;
    ++m_ndraw.m_local[icase];
}

void MagnitudeLogger::update_tau(double &tau, const ham_comp_t &gamma_sum) {
    if (m_static_tau) return;
    tau = m_max_bloom / gamma_sum;
    if (tau < m_tau_min) tau = m_tau_min;
    if (tau > m_tau_max) tau = m_tau_max;
}

void MagnitudeLogger::update(uint_t icycle, double &tau) {
    if (!icycle || icycle%m_period) return;
    m_gamma.all_max();
    update_tau(tau, m_gamma.m_reduced.sum());
}

void MagnitudeLogger::update(uint_t icycle, double &tau, ExcitGenGroup &excit_gens) {
    if (!icycle || icycle%m_period) return;
    m_gamma.all_max();
    auto gamma_sum = m_gamma.m_reduced.sum();
    // don't update any probabilities or timestep if all gammas are zero
    if (fptol::near_zero(gamma_sum)) return;
    update_tau(tau, gamma_sum);
    if (!m_static_probs){
        m_new_probs.assign(m_gamma.m_reduced.dbegin(), m_gamma.m_reduced.dend());
        for (auto& new_prob : m_new_probs) new_prob /= gamma_sum;
        prob::rectify(m_new_probs, m_prob_min);
        excit_gens.set_probs(m_new_probs);
        m_gamma.m_local.zero();
    }
}
