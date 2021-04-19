//
// Created by Robert John Anderson on 2020-04-11.
//

#include "StochasticPropagator.h"

void StochasticPropagator::add_boson_excitgen(const Hamiltonian<0> &ham) {}


void StochasticPropagator::add_boson_excitgen(const Hamiltonian<1> &ham) {
    m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
            new BosonExcitationGenerator(&ham, m_prng, ham.nboson_cutoff())));
}

StochasticPropagator::StochasticPropagator(const Hamiltonian<> &ham, const Options &opts, size_t npart) :
        Propagator(opts, ham, npart), m_prng(opts.prng_seed, opts.prng_ngen),
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


void StochasticPropagator::off_diagonal(Wavefunction &wf, const size_t &ipart) {
    const auto &row = wf.m_store.m_row;
    const defs::wf_t &weight = row.m_weight[ipart];
    ASSERT(!consts::float_is_zero(weight));
    ASSERT(consts::imag(weight) == 0.0 || m_ham.complex_valued())
    const auto &src_onv = row.m_onv;
    bool flag_initiator = row.m_initiator.get(ipart);
    bool flag_deterministic = row.m_deterministic.get(ipart);

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
        wf.add_spawn(m_dst_onv, delta, flag_initiator, flag_deterministic, ipart, src_onv, weight);
    }
}

void StochasticPropagator::diagonal(Wavefunction &wf, const size_t &ipart) {
    auto &row = wf.m_store.m_row;
    bool flag_deterministic = row.m_deterministic.get(ipart);
    const defs::ham_comp_t &hdiag = row.m_hdiag;
    if (flag_deterministic) {
        wf.scale_weight(ipart, 1 - (hdiag - m_shift[ipart]) * tau());
    } else {
        // the probability that each unit walker will die
        auto death_rate = (hdiag - m_shift[ipart]) * tau();
        if (death_rate <= 0.0 || death_rate > 1.0) {
            // clone  / create antiparticles continuously
            wf.scale_weight(ipart, 1 - death_rate);
        } else {
            // kill stochastically
            wf.set_weight(ipart, m_prng.stochastic_round(row.m_weight[ipart] * (1 - death_rate), m_opts.min_death_mag));
        }
    }
}
