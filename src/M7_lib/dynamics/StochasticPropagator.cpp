//
// Created by Robert John Anderson on 2020-04-11.
//

#include "StochasticPropagator.h"

StochasticPropagator::StochasticPropagator(const Hamiltonian &ham, const conf::Document &opts,
                                           const Wavefunction &wf) :
        Propagator(opts, ham, wf), m_prng(opts.m_prng.m_seed, opts.m_prng.m_ngen_block),
        m_excit_gen_group(ham, opts.m_propagator, m_prng),
        m_mag_log(opts.m_propagator.m_max_bloom,
                  opts.m_propagator.m_ndraw_min_for_dynamic,
                  m_excit_gen_group.ncase(),
                  opts.m_propagator.m_static_tau,
                  opts.m_propagator.m_static_probs,
                  opts.m_propagator.m_tau_min,
                  opts.m_propagator.m_tau_max,
                  opts.m_propagator.m_min_exlvl_prob,
                  opts.m_propagator.m_period),
        m_min_spawn_mag(opts.m_propagator.m_min_spawn_mag),
        m_min_death_mag(opts.m_propagator.m_min_death_mag){}


void StochasticPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    const defs::wf_t &weight = row.m_weight[ipart];
    /*
     * for bilinear estimators based on the consolidated annihilation of spawned contributions
     */
    defs::prob_t p_succeed_at_least_once = 1.0;

    DEBUG_ASSERT_NE(weight, 0.0, "should not attempt off-diagonal propagation from zero weight");
    DEBUG_ASSERT_TRUE(consts::imag(weight) == 0.0 || m_ham.complex_valued(),
                      "real-valued hamiltonian should never result in non-zero imaginary walker component")
    const auto &src_mbf = row.m_mbf;
    bool flag_initiator = row.m_initiator.get(ipart);
    bool flag_deterministic = row.m_deterministic.get(wf.iroot_part(ipart));

    const defs::wf_comp_t abs_weight = std::abs(weight);
    defs::prob_t prob_nattempt_floor, prob_gen, prob_thresh_accept;
    defs::ham_t helem;

    auto nattempt = size_t(m_prng.stochastic_round(abs_weight, 1.0, prob_nattempt_floor));
    if (!nattempt) return;
    /*
     * clear the cached-orbital representations of the current MBF
     */
    src_mbf.m_decoded.clear();

    auto &conn = m_conn[src_mbf];
    auto &dst_mbf = m_dst[src_mbf];
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {

        conn.clear();
        auto icase = m_excit_gen_group.draw_icase();
        if (!m_excit_gen_group.draw(icase, src_mbf, prob_gen, helem, conn)) {
            // null excitation generated
            continue;
        }
        m_mag_log.log(icase, helem, prob_gen);
        prob_gen *= m_excit_gen_group.get_prob(icase);

        conn.apply(src_mbf, dst_mbf);
        auto delta = -tau() * phase(weight) * helem / prob_gen;
        if (consts::nearly_zero(delta, 1e-14)) continue;
        imp_samp_delta(delta, src_mbf, dst_mbf, row.m_hdiag);
        /*
         * the stochastically-realized spawned contribution is equal to delta if delta is not lower in magnitude than
         * the minimum magnitude, otherwise it is stochastically rounded with respect to that magnitude
         */
        auto thresh_delta = m_prng.stochastic_threshold(delta, m_min_spawn_mag, prob_thresh_accept);
        if (thresh_delta==0.0) continue;

        if (wf.recv().m_row.m_send_parents) {
            /*
             * we require the parent weights to be sent along with the spawned contributions, so the instantaneous
             * weight must be "unbiased" by the probability that this connection was sampled a non-zero number of times.
             * there are 3 stochastic processes involved:
             *  1. number of attempts is either the *floor* or ceiling of the instantaneous weight (prob_nattempt_floor)
             *  2. generate a connected MBF (prob_gen)
             *  3. *accept* or reject based on the magnitude of the candidate walker (prob_thresh_accept)
             */
            defs::prob_t p_fail_one_attempt, p_fail_all_attempts;
            /*
             * failure to generate this connection at one attempt could be due to either
             *  1. generating some other connection via the excitation generator
             */
            p_fail_one_attempt = (1.0 - prob_gen);
            /*
             * 2. generating this connection only for it to be rejected by the stochastic threshold operation
             */
            p_fail_one_attempt += prob_gen * (1.0 - prob_thresh_accept);
            /*
             * now, failure to generate this connection at all attempts could be due to either
             *  1. deciding to draw the minimum number of attempts, only to fail on each attempt
             */
            auto floor = std::floor(abs_weight);
            p_fail_all_attempts = prob_nattempt_floor * std::pow(p_fail_one_attempt, floor);
            /*
             *  2. deciding to draw the maximum number of attempts, only to fail on each attempt
             */
            p_fail_all_attempts += (1.0 - prob_nattempt_floor) * std::pow(p_fail_one_attempt, floor+1.0);

            p_succeed_at_least_once = 1.0 - p_fail_all_attempts;
        }
        wf.add_spawn(dst_mbf, thresh_delta, flag_initiator, flag_deterministic,
                     ipart, src_mbf, weight / p_succeed_at_least_once);
    }
}

void StochasticPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    bool flag_deterministic = row.m_deterministic.get(wf.iroot_part(ipart));
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    if (flag_deterministic) {
        wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
    } else {
        // the probability that each unit walker will die
        auto death_rate = (hdiag - m_shift[ipart]) * tau();
        if (death_rate == 0.0) return;
        if (death_rate < 0.0 || death_rate > 1.0 || m_min_spawn_mag == 0.0) {
            // clone / create antiwalkers continuously
            wf.scale_weight(ipart, 1 - death_rate);
        } else {
            // kill stochastically
            wf.set_weight(ipart, m_prng.stochastic_round(row.m_weight[ipart] * (1 - death_rate), m_min_death_mag));
        }
    }
}

size_t StochasticPropagator::ncase_excit_gen() const {
    return m_excit_gen_group.ncase();
}

std::vector<defs::prob_t> StochasticPropagator::excit_gen_case_probs() const {
    return m_excit_gen_group.get_probs();
}

void StochasticPropagator::update(const size_t &icycle, const Wavefunction &wf) {
    Propagator::update(icycle, wf);
    m_mag_log.update(icycle, m_tau, m_excit_gen_group);
}
