//
// Created by rja on 04/08/2021.
//

#include "ExcitGenGroup.h"
#include "UniformFrmBos.h"

void ExcitGenGroup::init() {
    defs::prob_t norm = 0.0;
    for (size_t exsig=0ul; exsig<defs::nexsig; ++exsig){
        const auto ptr = m_exgens[exsig].get();
        if (!ptr) continue;
        if (conn_utils::pure_frm(exsig)) m_frm_inds.push_back(m_active_exsigs.size());
        m_probs.push_back(ptr->approx_nconn());
        norm+=m_probs.back();
        m_active_exsigs.push_back(exsig);
    }
    for (auto &prob: m_probs) prob /= norm;
    update_cumprobs();
}

void ExcitGenGroup::update_cumprobs() {
    m_cumprobs = m_probs;
    m_frm_probs.clear();
    m_frm_norm = 0.0;
    auto it = m_probs.cbegin();
    for (const auto& exsig: m_active_exsigs){
        if (conn_utils::pure_frm(exsig)) {
            m_frm_probs.push_back(*it);
            m_frm_norm+=m_frm_probs.back();
        }
        ++it;
    }

    for (size_t iprob = 1ul; iprob < m_probs.size(); ++iprob)
        m_cumprobs[iprob] = m_cumprobs[iprob - 1] + m_probs[iprob];
    DEBUG_ASSERT_TRUE(consts::floats_nearly_equal(1.0, m_cumprobs.back()),
                      "cumulative probability should be 1.0");
    m_frm_cumprobs = m_frm_probs;
    for (size_t iprob = 1ul; iprob < m_frm_probs.size(); ++iprob)
        m_frm_cumprobs[iprob] = m_frm_cumprobs[iprob - 1] + m_frm_probs[iprob];
    DEBUG_ASSERT_TRUE(consts::floats_nearly_equal(1.0, m_frm_cumprobs.back()),
                      "cumulative probability should be 1.0");
}

ExcitGenGroup::ExcitGenGroup(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng) :
        m_prng(prng) {
    if (ham.m_frm.is_hubbard_1d() || ham.m_frm.is_hubbard_1d_pbc()) {
        m_exgens[conn_utils::exsig(1,1,0,0)] = std::unique_ptr<ExcitGen>(
                new Hubbard1dSingles(ham, prng, ham.m_frm.is_hubbard_1d_pbc()));
    } else {
        m_exgens[conn_utils::exsig(1,1,0,0)] = std::unique_ptr<ExcitGen>(
                new UniformSingles(ham, prng));
    }
    if (ham.m_frm.int_2e_rank()) {
        if (opts.m_excit_gen.get() == "pchb") {
            m_exgens[conn_utils::exsig(2,2,0,0)] = std::unique_ptr<ExcitGen>(
                    new HeatBathDoubles(ham, prng));
        }
    }
    if (ham.m_bos.m_nboson_max) {
        m_exgens[conn_utils::exsig(0,0,1,0)] =
                std::unique_ptr<ExcitGen>(new UniformFrmBos(ham, prng, true));
        m_exgens[conn_utils::exsig(0,0,0,1)] =
                std::unique_ptr<ExcitGen>(new UniformFrmBos(ham, prng, false));
    }

    init();

    auto spec_probs = opts.m_exlvl_probs_init.get();
    if (spec_probs.size()) {
        log::info("Using probabilities specified in config to initialize excitation generator group");
        REQUIRE_TRUE_ALL(spec_probs.empty() || spec_probs.size() == size() ||
                         spec_probs.size() == size() - 1, "invalid number of probabilities specified in config");
        if (spec_probs.size() + 1 == size()) {
            auto norm = std::accumulate(spec_probs.cbegin(), spec_probs.cend(), 0.0);
            REQUIRE_LE_ALL(norm, 1.0, "no valid value of the remaining probability can normalize the distribution");
            spec_probs.push_back(1.0 - norm);
        }
        set_probs(spec_probs);
    }
}

size_t ExcitGenGroup::size() const {
    return m_active_exsigs.size();
}

void ExcitGenGroup::set_probs(const std::vector<defs::prob_t> &probs) {
    DEBUG_ASSERT_EQ(probs.size(), size(), "incorrect number of probabilities given");
    m_probs = probs;
    update_cumprobs();
}

defs::prob_t ExcitGenGroup::get_prob(const size_t &iex) const {
    DEBUG_ASSERT_LT(iex, size(), "excit gen index OOB");
    return m_probs[iex];
}

defs::prob_t ExcitGenGroup::get_prob_frm(const size_t &iex) const {
    return get_prob(iex) / m_frm_norm;
}

const std::vector<defs::prob_t> &ExcitGenGroup::get_probs() const {
    return m_probs;
}

ExcitGen &ExcitGenGroup::operator[](const size_t &iex) {
    DEBUG_ASSERT_LT(iex, size(), "excit gen index OOB");
    return *m_exgens[m_active_exsigs[iex]];
}

const ExcitGen &ExcitGenGroup::operator[](const size_t &iex) const {
    DEBUG_ASSERT_LT(iex, size(), "excit gen index OOB");
    return *m_exgens[m_active_exsigs[iex]];
}

size_t ExcitGenGroup::draw_iex() {
    auto r = m_prng.draw_float();
    for (size_t i = 0ul; i < size(); ++i) {
        if (r < m_cumprobs[i]) return i;
    }
    return size() - 1;
}

size_t ExcitGenGroup::draw_iex_frm() {
    auto r = m_prng.draw_float();
    auto size = m_frm_inds.size();
    for (size_t i = 0ul; i < size; ++i) {
        if (r < m_frm_cumprobs[i]) return m_frm_inds[i];
    }
    return m_frm_inds[size - 1];
}

void ExcitGenGroup::log_breakdown() const {
    log::info("Excitation class probability breakdown:");
    for (size_t i = 0ul; i < size(); ++i)
        log::info("{:<40} {}", (*this)[i].description(), get_prob(i));
}
