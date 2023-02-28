//
// Created by Robert John Anderson on 2020-04-11.
//

#include "StochLinear.h"

StochLinear::StochLinear(const Hamiltonian& ham, const conf::Document& opts,
                                           const wf::Fci& wf) :
        Propagator(opts, ham, wf), m_prng(opts.m_prng.m_seed, opts.m_prng.m_ngen_block),
        m_excit_gen_group(ham, opts.m_propagator, m_prng, wf.m_sector.particles()),
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
        m_min_death_mag(opts.m_propagator.m_min_death_mag){
}


void StochLinear::diagonal(wf::Fci& wf, Walker& walker, const uint_t& ipart) {
    bool flag_deterministic = walker.m_deterministic.get(wf.iroot_part(ipart));
    const ham_comp_t& hdiag = walker.m_hdiag;
    if (flag_deterministic) {
        wf.scale_weight(walker, ipart, 1 - (hdiag - m_shift[ipart]) * tau());
    } else {
        // the probability that each unit walker will die
        auto death_rate = (hdiag - m_shift[ipart]) * tau();
        if (death_rate == 0.0) return;
        if (death_rate < 0.0 || death_rate > 1.0 || m_min_death_mag == 0.0) {
            // clone / create antiwalkers continuously
            wf.scale_weight(walker, ipart, 1 - death_rate);
        } else {
            // kill stochastically
            auto new_weight = m_prng.stochastic_round(walker.m_weight[ipart] * (1 - death_rate), m_min_death_mag);
            wf.set_weight(walker, ipart, new_weight);
        }
    }
}

void StochLinear::off_diagonal(wf::Fci& wf, const Walker& walker, const uint_t& ipart) {
    const wf_t& weight = walker.m_weight[ipart];
    /*
     * for bilinear estimators based on the consolidated annihilation of spawned contributions
     */
    prob_t p_succeed_at_least_once = 1.0;

    DEBUG_ASSERT_NE(weight, 0.0, "should not attempt off-diagonal propagation from zero weight");
    DEBUG_ASSERT_TRUE(m_ham.complex_valued() || fptol::near_real(weight),
                      "real-valued hamiltonian should never result in non-zero imaginary walker component")
    const auto& src_mbf = walker.m_mbf;
    bool is_initiator = walker.is_initiator(ipart, m_nadd_initiator);
    bool flag_deterministic = walker.m_deterministic.get(wf.iroot_part(ipart));

    const wf_comp_t abs_weight = std::abs(weight);
    prob_t prob_nattempt_ceil, prob_gen, prob_thresh_accept;
    ham_t helem;

    auto nattempt = uint_t(m_prng.stochastic_round(abs_weight, 1.0, prob_nattempt_ceil));
    if (!nattempt) return;
    /*
     * clear the cached-orbital representations of the current MBF
     */
    src_mbf.m_decoded.clear();

    auto& conn = m_conn[src_mbf];
    auto& dst_mbf = m_dst[src_mbf];
    for (uint_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {

        conn.clear();
        auto icase = m_excit_gen_group.draw_icase();
        if (!m_excit_gen_group.draw(icase, src_mbf, prob_gen, helem, conn)) {
            // null excitation generated
            continue;
        }
        /*
         * the magnitude logging is case-resolved, so the relevant magnitude is the one scaled by the reciprocal of the
         * probability that conn was drawn from src given that the case was picked
         */
        m_mag_log.log(icase, helem, prob_gen);
        /*
         * scale by the probability of drawing this case (and add probs of other cases which could have produced the
         * same connection if required)
         */
        m_excit_gen_group.update_prob(icase, src_mbf, prob_gen, conn);

        conn.apply(src_mbf, dst_mbf);
        auto delta = -tau() * phase(weight) * helem / prob_gen;
        if (fptol::near_zero(delta)) continue;
        imp_samp_delta(delta, src_mbf, dst_mbf);
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
             *  1. number of attempts is either the floor or ceiling of the instantaneous weight (prob_nattempt_ceil)
             *  2. generate a connected MBF (prob_gen)
             *  3. *accept* or reject based on the magnitude of the candidate walker (prob_thresh_accept)
             */
            prob_t p_fail_one_attempt, p_fail_all_attempts;
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
             *  1. deciding to draw the maximum number of attempts, only to fail on each attempt
             */
            auto ceil = std::ceil(abs_weight);
            p_fail_all_attempts = prob_nattempt_ceil * std::pow(p_fail_one_attempt, ceil);
            /*
             *  2. deciding to draw the minimum number of attempts, only to fail on each attempt
             */
            p_fail_all_attempts += (1.0 - prob_nattempt_ceil) * std::pow(p_fail_one_attempt, ceil - 1.0);

            p_succeed_at_least_once = 1.0 - p_fail_all_attempts;
        }
        wf.add_spawn(dst_mbf, thresh_delta, is_initiator, flag_deterministic,
                     ipart, src_mbf, weight / p_succeed_at_least_once);
    }
}

uint_t StochLinear::ncase_excit_gen() const {
    return m_excit_gen_group.ncase();
}

v_t<prob_t> StochLinear::excit_gen_case_probs() const {
    return m_excit_gen_group.get_probs();
}

void StochLinear::update(uint_t icycle, const wf::Fci& wf, const wf::Refs& refs) {
    Propagator::update(icycle, wf, refs);
    m_mag_log.update(icycle, m_tau, m_excit_gen_group);
}

hash::digest_t StochLinear::checksum_() const {
    return m_prng.checksum();
}
