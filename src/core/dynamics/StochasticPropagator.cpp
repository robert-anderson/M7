//
// Created by Robert John Anderson on 2020-04-11.
//

#if 0
#include "StochasticPropagator.h"

void StochasticPropagator::add_boson_excitgen(const Hamiltonian<0> &ham) {}


void StochasticPropagator::add_boson_excitgen(const Hamiltonian<1> &ham) {
    m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
            new BosonExcitationGenerator(&ham, m_prng, ham.nboson_cutoff())));
}

StochasticPropagator::StochasticPropagator(const Hamiltonian<> &ham, const Options &opts) :
        Propagator(opts, ham), m_prng(opts.prng_seed, opts.prng_ngen),
        m_min_spawn_mag(opts.min_spawn_mag) {

    m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
            new UniformSingles(&m_ham, m_prng)));
    if (ham.int_2e_rank() && opts.excit_gen == "pchb") {
        m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
                new HeatBathDoubles(&m_ham, m_prng)));
    }
    add_boson_excitgen(ham);

    m_exgen_drawer = std::unique_ptr<WeightedDrawer>(new WeightedDrawer(m_exgens.size(), m_prng));

    const defs::prob_t prob_boson = 0.2; // TODO: make dynamic.
    if (m_exgens.size() == 2)
        m_exgen_drawer->set(m_magnitude_logger.m_psingle);
    else if (m_exgens.size() == 3)
        m_exgen_drawer->set(m_magnitude_logger.m_psingle, 1.0 - m_magnitude_logger.m_psingle - prob_boson);

    log::info("Excitation class probability breakdown {}", utils::to_string(m_exgen_drawer->m_probs));
}


void StochasticPropagator::off_diagonal(Wavefunction &m_wf, const size_t &irow) {
    auto weight = m_wf.m_store.m_weight(irow, 0, 0);
    ASSERT(!consts::float_is_zero(weight));
    ASSERT(consts::imag(weight) == 0.0 || m_ham.complex_valued())
    auto src_onv = m_wf.m_store.m_onv(irow);
    bool flag_initiator = m_wf.m_store.m_flags.m_initiator(irow, 0, 0);
    bool flag_deterministic = m_wf.m_store.m_flags.m_deterministic(irow);

    m_occ.update(src_onv);
    m_vac.update(src_onv);
    size_t nattempt = get_nattempt(weight);
#ifdef VERBOSE_DEBUGGING
    std::cout << consts::verb << "spawn attempts: " << nattempt << std::endl;
#endif
    defs::prob_t prob;
    defs::ham_t helem;
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        m_aconn.zero();
        m_dst_onv = src_onv;
        size_t iexgen = m_exgen_drawer->draw();
        auto &exgen = m_exgens[iexgen];
        if (!exgen->draw(src_onv, m_dst_onv, m_occ, m_vac, prob, helem, m_aconn)) continue;
        prob *= m_exgen_drawer->prob(iexgen);
        auto delta = -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob;
        if (consts::float_is_zero(delta)) continue;
        delta = m_prng.stochastic_threshold(delta, m_opts.min_spawn_mag);
        if (consts::float_is_zero(delta)) continue;
        m_wf.add_spawn(m_dst_onv, delta, flag_initiator, flag_deterministic);
    }
}

void StochasticPropagator::diagonal(Wavefunction &m_wf, const size_t &irow) {
    bool flag_deterministic = m_wf.m_store.m_flags.m_deterministic(irow);
    auto hdiag = m_wf.m_store.m_hdiag(irow);
    if (flag_deterministic) {
        m_wf.scale_weight(irow, 1 - (hdiag - m_shift) * tau());
    } else {
        // the probability that each unit walker will die
        auto death_rate = (hdiag - m_shift) * tau();
        if (death_rate < 0.0 || death_rate > 1.0) {
            // clone  / create antiparticles continuously
            m_wf.scale_weight(irow, 1 - death_rate);
        } else {
            auto weight = m_wf.m_store.m_weight(irow, 0, 0);
            // kill stochastically
            m_wf.set_weight(irow, m_prng.stochastic_round(weight * (1 - death_rate), m_opts.min_death_mag));
        }
    }
}
#endif