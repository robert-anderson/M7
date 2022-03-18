//
// Created by Robert John Anderson on 2020-04-11.
//

#include <excitgen/HubbardUniform.h>
#include "StochasticPropagator.h"


StochasticPropagator::StochasticPropagator(const Hamiltonian &ham, const fciqmc_config::Document &opts,
                                           const NdFormat<defs::ndim_wf> &wf_fmt) :
        Propagator(opts, ham, wf_fmt), m_prng(opts.m_prng.m_seed, opts.m_prng.m_ngen_block),
        m_excit_gens(ham, opts.m_propagator, m_prng),
        m_mag_log(opts.m_propagator.m_max_bloom,
                  opts.m_propagator.m_ndraw_min_for_dynamic,
                  m_excit_gens.size(),
                  opts.m_propagator.m_static_tau,
                  opts.m_propagator.m_static_probs,
                  opts.m_propagator.m_tau_min,
                  opts.m_propagator.m_tau_max,
                  opts.m_propagator.m_min_exlvl_prob,
                  opts.m_propagator.m_period),
        m_min_spawn_mag(opts.m_propagator.m_min_spawn_mag) {
    m_excit_gens.log_breakdown();
}


void StochasticPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    const defs::wf_t &weight = row.m_weight[ipart];
    double rdm_factor = 1.0;

    DEBUG_ASSERT_NE(weight, 0.0, "should not attempt offdiagonal propagation from zero weight");
    DEBUG_ASSERT_TRUE(consts::imag(weight) == 0.0 || m_ham.complex_valued(),
                      "real-valued hamiltonian should never result in non-zero imaginary walker component")
    const auto &src_mbf = row.m_mbf;
    bool flag_initiator = row.m_initiator.get(ipart);
    bool flag_deterministic = row.m_deterministic.get(wf.iroot_part(ipart));

    size_t nattempt = get_nattempt(weight);
    if (!nattempt) return;

    m_excit_gens.clear_cached_orbs();

    defs::prob_t prob;
    defs::ham_t helem;

    auto &conn = m_conn[src_mbf];
    auto &dst_mbf = m_dst[src_mbf];
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {

        conn.clear();
        auto iex = m_excit_gens.draw_iex();
        if (!m_excit_gens.draw(iex, src_mbf, prob, helem, conn)) {
            // null excitation generated
            continue;
        }
        m_mag_log.log(iex, helem, prob);
        prob *= m_excit_gens.get_prob(iex);

        conn.apply(src_mbf, dst_mbf);
        auto delta = -tau() * phase(weight) * helem / prob;
        if (consts::nearly_zero(delta, 1e-14)) continue;
        imp_samp_delta(delta, src_mbf, dst_mbf, row.m_hdiag);
        delta = m_prng.stochastic_threshold(delta, m_opts.m_propagator.m_min_spawn_mag);
        if (consts::nearly_zero(delta, 1e-14)) continue;

        if (wf.recv().m_row.m_send_parents) {
            // reweight by probability that this connection was sampled a non-zero number of times
            rdm_factor = 1.0 / (1.0 - std::pow(1 - prob, nattempt));
        }
        wf.add_spawn(dst_mbf, delta, flag_initiator, flag_deterministic,
                     ipart, src_mbf, rdm_factor * weight);
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
        if (death_rate < 0.0 || death_rate > 1.0 || m_opts.m_propagator.m_min_spawn_mag == 0.0) {
            // clone / create antiwalkers continuously
            wf.scale_weight(ipart, 1 - death_rate);
        } else {
            // kill stochastically
            wf.set_weight(ipart, m_prng.stochastic_round(row.m_weight[ipart] * (1 - death_rate),
                                                         m_opts.m_propagator.m_min_death_mag));
        }
    }
}

size_t StochasticPropagator::nexcit_gen() const {
    return m_excit_gens.size();
}

std::vector<defs::prob_t> StochasticPropagator::exlvl_probs() const {
    return m_excit_gens.get_probs();
}

void StochasticPropagator::update(const size_t &icycle, const Wavefunction &wf) {
    Propagator::update(icycle, wf);
    m_mag_log.update(icycle, m_tau, m_excit_gens);
}
