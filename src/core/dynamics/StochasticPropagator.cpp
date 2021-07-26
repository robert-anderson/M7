//
// Created by Robert John Anderson on 2020-04-11.
//

#include <src/core/excitgen/HubbardSingles.h>
#include "StochasticPropagator.h"

void StochasticPropagator::add_boson_excitgen(const Hamiltonian<0> &ham) {}


void StochasticPropagator::add_boson_excitgen(const Hamiltonian<1> &ham) {
//    m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
//            new BosonExcitationGenerator(&ham, m_prng, ham.nboson_cutoff())));
}

StochasticPropagator::StochasticPropagator(const Hamiltonian<> &ham, const fciqmc_config::Document &opts, const NdFormat<defs::ndim_wf>& wf_fmt):
        Propagator(opts, ham, wf_fmt), m_prng(opts.m_prng.m_seed, opts.m_prng.m_ngen_block),
        m_min_spawn_mag(opts.m_propagator.m_min_spawn_mag) {


    if (ham.is_hubbard() || ham.is_hubbard_pbc()){
        m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
                new HubbardSingles(&m_ham, m_prng, ham.is_hubbard_pbc())));
    }
    else {
        m_exgens.push_back(std::unique_ptr<ExcitationGenerator>(
                new UniformSingles(&m_ham, m_prng)));
    }
    if (ham.int_2e_rank() && opts.m_propagator.m_excit_gen.get() == "pchb") {
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
    const defs::wf_t& weight = row.m_weight[ipart];
    double rdm_unbias_factor = 1.0;

    ASSERT(!consts::float_is_zero(weight));
    ASSERT(consts::imag(weight) == 0.0 || m_ham.complex_valued())
    const auto &src_onv = row.m_mbf;
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

    auto& conn = m_conn[src_onv];
    auto& dst_onv = m_dst[src_onv];
    for (size_t iattempt = 0ul; iattempt < nattempt; ++iattempt) {
        conn.clear();
        dst_onv = src_onv;
        size_t iexgen = m_exgen_drawer->draw();
        auto &exgen = m_exgens[iexgen];
        if (!exgen->draw(src_onv, dst_onv, m_occ, m_vac, prob, helem, conn)) continue;
        prob *= m_exgen_drawer->prob(iexgen);
        auto delta = -(weight / (defs::ham_comp_t) nattempt) * tau() * helem / prob;
        if (consts::float_is_zero(delta)) continue;
        delta = m_prng.stochastic_threshold(delta, m_opts.m_propagator.m_min_spawn_mag);
        if (consts::float_is_zero(delta)) continue;

        if (wf.recv().m_row.m_send_parents){
            if (m_opts.m_propagator.m_consolidate_spawns) {
                // reweight by probability that this connection was sampled a non-zero number of times
                rdm_unbias_factor = 1.0 / (1.0 - std::pow(1 - prob, nattempt));
            }
            else {
                // reweight for expected number of draws of this connection
                rdm_unbias_factor = 1.0 / (prob*nattempt);
            }
        }
        wf.add_spawn(dst_onv, delta, flag_initiator, flag_deterministic,
                     ipart, src_onv, rdm_unbias_factor*weight);
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
        if (death_rate==0.0) return;
        if (death_rate < 0.0 || death_rate > 1.0 || m_opts.m_propagator.m_min_spawn_mag==0.0) {
            // clone / create antiwalkers continuously
            wf.scale_weight(ipart, 1 - death_rate);
        } else {
            // kill stochastically
            wf.set_weight(ipart, m_prng.stochastic_round(row.m_weight[ipart] * (1 - death_rate),
                                                         m_opts.m_propagator.m_min_death_mag));
        }
    }
}

bool StochasticPropagator::is_exact() const {
    return false;
}
