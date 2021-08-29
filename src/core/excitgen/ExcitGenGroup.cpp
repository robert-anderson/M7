//
// Created by rja on 04/08/2021.
//

#include "ExcitGenGroup.h"
#include "LadderPureHolstein.h"

void ExcitGenGroup::init() {
    defs::prob_t norm = 0.0;
    for (size_t exsig = 0ul; exsig < defs::nexsig; ++exsig) {
        auto ptr = m_exsigs_to_exgens[exsig];
        if (!ptr) continue;
        if (exsig_utils::is_pure_frm(exsig)) m_frm_inds.push_back(m_active_exsigs.size());
        m_probs.push_back(ptr->approx_nconn());
        norm += m_probs.back();
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
    for (const auto &exsig: m_active_exsigs) {
        if (exsig_utils::is_pure_frm(exsig)) {
            m_frm_probs.push_back(*it);
            m_frm_norm += m_frm_probs.back();
        }
        ++it;
    }
    for (auto& frm_prob: m_frm_probs) frm_prob/=m_frm_norm;

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

void ExcitGenGroup::add(std::unique_ptr<ExcitGen> &&exgen, const inds &exsigs) {
    m_exgens.emplace_front(std::move(exgen));
    auto ptr = m_exgens.begin()->get();
    for (auto &exsig: exsigs) {
        REQUIRE_LT(exsig, defs::nexsig, "exsig OOB");
        REQUIRE_TRUE(m_exsigs_to_exgens[exsig] == nullptr,
                     "can't specify more than one excitation generator for the same exsig in an ExcitGenGroup");
        m_exsigs_to_exgens[exsig] = ptr;
    }
}

void ExcitGenGroup::add(std::unique_ptr<ExcitGen> &&exgen) {
    auto exsigs = exgen->m_exsigs;
    add(std::move(exgen), exsigs);
}

ExcitGenGroup::ExcitGenGroup(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng) :
        m_prng(prng), m_cached_orbs(ham.m_frm.m_point_group_map) {
    m_exsigs_to_exgens.fill(nullptr);
    bool any_singles =
            ham.m_frm.m_contribs_1100.is_nonzero(ex_single) || ham.m_frm.m_contribs_2200.is_nonzero(ex_single);
    if (any_singles) {
        if (ham.m_frm.m_model_attrs.is_hubbard_1d() || ham.m_frm.m_model_attrs.is_hubbard_1d_pbc()) {
            add(std::unique_ptr<ExcitGen>(new Hubbard1dSingles(ham, prng)));
        } else {
            add(std::unique_ptr<ExcitGen>(new UniformSingles(ham, prng)));
        }
    }
    if (ham.m_frm.m_contribs_2200.is_nonzero(ex_double)) {
        if (opts.m_excit_gen.get() == "pchb")
            add(std::unique_ptr<ExcitGen>(new HeatBathDoubles(ham, prng)));
    }
    if (ham.m_bos.m_nboson_max) {
        add(std::unique_ptr<ExcitGen>(new LadderPureUniform(ham, prng)),
            {exsig_utils::ex_0010, exsig_utils::ex_0001});
//        m_exgens[conn_utils::encode_exsig(0, 0, 1, 0)] =
//                std::unique_ptr<ExcitGen>(new UniformHolstein(ham, prng, true));
//        m_exgens[conn_utils::encode_exsig(0, 0, 0, 1)] =
//                std::unique_ptr<ExcitGen>(new UniformHolstein(ham, prng, false));
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
    return *m_exsigs_to_exgens[m_active_exsigs[iex]];
}

const ExcitGen &ExcitGenGroup::operator[](const size_t &iex) const {
    DEBUG_ASSERT_LT(iex, size(), "excit gen index OOB");
    return *m_exsigs_to_exgens[m_active_exsigs[iex]];
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
    log::info("Excitation generation probability breakdown by excitation signature:");
    for (size_t i = 0ul; i < size(); ++i) {
        auto exsig = exsig_utils::to_string(m_active_exsigs[i]);
        log::info("excitation generator {:<40} for exsig {}:  prob={}",
                  (*this)[i].description(), exsig, get_prob(i));
    }
}

bool ExcitGenGroup::draw(const size_t &iex, const FrmOnv &src, prob_t &prob, ham_t &helem, conn::FrmOnv &conn) {
    auto exsig = m_active_exsigs[iex];
    return (*this)[iex].draw(exsig, src, m_cached_orbs, prob, helem, conn);
}

bool ExcitGenGroup::draw(const size_t &iex, const FrmBosOnv &src, prob_t &prob, ham_t &helem, conn::FrmBosOnv &conn) {
    auto exsig = m_active_exsigs[iex];
    return (*this)[iex].draw(exsig, src, m_cached_orbs, prob, helem, conn);
}

void ExcitGenGroup::clear_cached_orbs() {
    m_cached_orbs.clear();
}