//
// Created by rja on 06/04/2022.
//

#include "ExcitGenGroup2.h"

void ExcitGenGroup2::update_cumprobs() {
    m_cumprobs = m_probs;
    for (size_t icase = 1ul; icase < m_probs.size(); ++icase)
        m_cumprobs[icase] = m_cumprobs[icase - 1] + m_probs[icase];
    DEBUG_ASSERT_NEARLY_EQ(m_cumprobs.back(), 1.0, consts::eps<prob_t>(), "cumulative probability should be 1.0");
}

ExcitGenGroup2::ExcitGenGroup2(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng) :
        m_prng(prng) {
    ExcitGen2::excit_gen_list_t list;
    list = ham.m_frm->make_excit_gens(prng, opts);
    m_list.merge(list);
    list = ham.m_frmbos->make_excit_gens(prng, opts);
    m_list.merge(list);
    list = ham.m_bos->make_excit_gens(prng, opts);
    m_list.merge(list);
    REQUIRE_FALSE(m_list.empty(), "there should be at least one active excitation case");
    /*
     * exsigs have a many-to-many relationship with excitation generators.
     */
    m_probs.clear();
    defs::prob_t norm = 0.0;
    for (const auto& excit_gen: m_list) {
        for (auto& exsig: excit_gen->m_exsigs){
            m_excit_cases.push_back({exsig, excit_gen.get()});
            m_probs.push_back(defs::prob_t(m_excit_cases.back().m_excit_gen->approx_nconn()));
            norm+=m_probs.back();
        }
    }
    for (auto &prob: m_probs) prob /= norm;
    update_cumprobs();
}

size_t ExcitGenGroup2::ncase() const {
    return m_excit_cases.size();
}

size_t ExcitGenGroup2::draw_icase() {
    auto r = m_prng.draw_float();
    for (size_t i = 0ul; i < ncase(); ++i) {
        if (r < m_cumprobs[i]) return i;
    }
    return ncase() - 1;
}

ExcitCase &ExcitGenGroup2::operator[](size_t icase) {
    return m_excit_cases[icase];
}

void ExcitGenGroup2::set_probs(const std::vector<defs::prob_t> &probs) {
    DEBUG_ASSERT_EQ(probs.size(), ncase(), "incorrect number of probabilities given");
    m_probs = probs;
    update_cumprobs();
}

defs::prob_t ExcitGenGroup2::get_prob(const size_t &icase) const {
    DEBUG_ASSERT_LT(icase, ncase(), "excit gen case index OOB");
    return m_probs[icase];
}

const std::vector<defs::prob_t>& ExcitGenGroup2::get_probs() const {
    return m_probs;
}

bool ExcitGenGroup2::draw(const size_t &icase, const FrmOnv &src, prob_t &prob, ham_t &helem, conn::FrmOnv &conn) {
    auto& excase = m_excit_cases[icase];
    return excase.m_excit_gen->draw(excase.m_exsig, src, prob, helem, conn);
}

bool ExcitGenGroup2::draw(const size_t &icase, const FrmBosOnv &src, prob_t &prob, ham_t &helem, conn::FrmBosOnv &conn) {
    auto& excase = m_excit_cases[icase];
    return excase.m_excit_gen->draw(excase.m_exsig, src, prob, helem, conn);
}

bool ExcitGenGroup2::draw(const size_t &icase, const BosOnv &src, prob_t &prob, ham_t &helem, conn::BosOnv &conn) {
    auto& excase = m_excit_cases[icase];
    return excase.m_excit_gen->draw(excase.m_exsig, src, prob, helem, conn);
}

void ExcitGenGroup2::log() const {
    std::vector<std::vector<std::string>> rows = {{"Excitation Signature", "Description", "Probability"}};
    for (size_t icase=0ul; icase<ncase(); ++icase){
        auto exsig_str = exsig_utils::to_string(m_excit_cases[icase].m_exsig);
        auto prob_str = utils::to_string(get_prob(icase));
        rows.push_back({exsig_str, m_excit_cases[icase].m_excit_gen->m_description, prob_str});
    }
    log::info_table("Excitation generation breakdown", rows, true);
}
