//
// Created by rja on 04/08/2021.
//

#include "ExcitGenGroup.h"

void ExcitGenGroup::init_probs() {
    for (auto ptr: m_ptrs) m_probs.push_back(ptr->approx_nconn());
    auto norm = std::accumulate(m_probs.cbegin(), m_probs.cend(), 0.0);
    for (auto &prob: m_probs) prob /= norm;
    update_cumprobs();
}

void ExcitGenGroup::update_cumprobs() {
    m_cumprobs = m_probs;
    for (size_t iprob = 1ul; iprob < m_probs.size(); ++iprob)
        m_cumprobs[iprob] = m_cumprobs[iprob - 1] + m_probs[iprob];
    DEBUG_ASSERT_TRUE(consts::floats_nearly_equal(1.0, m_cumprobs.back()),
                      "cumulative probability should be 1.0");
}

ExcitGenGroup::ExcitGenGroup(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng) :
        m_prng(prng) {
    if (ham.m_frm.is_hubbard_1d() || ham.m_frm.is_hubbard_1d_pbc()) {
        m_frm_singles = std::unique_ptr<FrmExcitGen>(
                new Hubbard1dSingles(ham, prng, ham.m_frm.is_hubbard_1d_pbc()));
    } else {
        m_frm_singles = std::unique_ptr<FrmExcitGen>(
                new UniformSingles(ham, prng));
    }
    if (ham.m_frm.int_2e_rank()) {
        if (opts.m_excit_gen.get() == "pchb") {
            m_frm_doubles = std::unique_ptr<FrmExcitGen>(
                    new HeatBathDoubles(ham, prng));
        }
    }
    if (ham.m_bos.m_nboson_max) {
        m_frmbos = std::unique_ptr<FrmBosExcitGen>(
                new FrmBosExcitGen(ham, prng));
    }

    // add to ptrs in reverse order of generally expected precedence:
    if (m_frm_doubles) m_ptrs.push_back(m_frm_doubles.get());
    if (m_frm_singles) m_ptrs.push_back(m_frm_singles.get());
    if (m_frmbos) m_ptrs.push_back(m_frmbos.get());
    init_probs();

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
    return m_ptrs.size();
}

void ExcitGenGroup::set_probs(const std::vector<defs::prob_t> &probs) {
    DEBUG_ASSERT_EQ(probs.size(), size(), "incorrect number of probabilities given");
    m_probs = probs;
    update_cumprobs();
}

void ExcitGenGroup::set_prob(const defs::prob_t &prob, ExcitGen *ptr) {
    DEBUG_ASSERT_LE(prob, 1.0, "given prob exceeds one");
    DEBUG_ASSERT_GE(prob, 0.0, "given prob is less than zero");
    auto norm_remain = 1.0 - prob;
    for (size_t i = 0ul; i < size(); ++i) {
        if (m_ptrs[i] != ptr) m_probs[i] *= norm_remain;
    }
    update_cumprobs();
}

const defs::prob_t &ExcitGenGroup::get_prob(const size_t &iexlvl) const {
    DEBUG_ASSERT_LT(iexlvl, size(), "excit gen index OOB");
    return m_probs[iexlvl];
}

const std::vector<defs::prob_t> &ExcitGenGroup::get_probs() const {
    return m_probs;
}

ExcitGen &ExcitGenGroup::operator[](const size_t &iexlvl) {
    DEBUG_ASSERT_LT(iexlvl, size(), "excit gen index OOB");
    return *m_ptrs[iexlvl];
}

size_t ExcitGenGroup::draw_iexlvl() {
    auto r = m_prng.draw_float();
    for (size_t i = 0ul; i < size(); ++i) {
        if (r < m_cumprobs[i]) return i;
    }
    return size() - 1;
}

void ExcitGenGroup::log_breakdown() const {
    log::info("Excitation class probability breakdown:");
    for (size_t i = 0ul; i < size(); ++i)
        log::info("{:<40} {}", m_ptrs[i]->description(), get_prob(i));
}
