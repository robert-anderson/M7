//
// Created by Robert J. Anderson on 06/04/2022.
//

#include "ExcitGenGroup.h"

void ExcitGenGroup::update_cumprobs() {
    m_cumprobs = m_probs;
    for (uint_t icase = 1ul; icase < m_probs.size(); ++icase)
        m_cumprobs[icase] = m_cumprobs[icase - 1] + m_probs[icase];
    DEBUG_ASSERT_NUM_EQ(m_cumprobs.back(), 1.0, "cumulative probability should be 1.0");
}

ExcitGenGroup::ExcitGenGroup(const Hamiltonian& ham, const conf::Propagator& opts, PRNG& prng) :
        m_prng(prng){
    ExcitGen::excit_gen_list_t list;
    list = ham.m_frm.make_excit_gens(prng, opts);
    m_list.merge(list);
    list = ham.m_frmbos.make_excit_gens(prng, opts);
    m_list.merge(list);
    list = ham.m_bos.make_excit_gens(prng, opts);
    m_list.merge(list);
    REQUIRE_FALSE(m_list.empty(), "there should be at least one active excitation case");
    /*
     * exsigs have a many-to-many relationship with excitation generators.
     */
    m_probs.clear();

    for (const auto& excit_gen: m_list) {
        for (auto& exsig: excit_gen->m_exsigs){
            m_excit_cases.push_back({exsig, excit_gen.get()});
            // let all cases be equally likely
            m_probs.push_back(1.0);
        }
    }
    m_exsig_icases.resize(exsig::c_ndistinct, uintv_t());
    // fill the map from exsigs to exgens
    for (uint_t icase=0ul; icase<m_excit_cases.size(); ++icase) m_exsig_icases[m_excit_cases[icase].m_exsig].push_back(icase);
    set_probs(m_probs);
    log();
}

uint_t ExcitGenGroup::ncase() const {
    return m_excit_cases.size();
}

uint_t ExcitGenGroup::draw_icase() {
    auto r = m_prng.draw_float();
    for (uint_t i = 0ul; i < ncase(); ++i) {
        if (r < m_cumprobs[i]) return i;
    }
    return ncase() - 1;
}

ExcitCase& ExcitGenGroup::operator[](uint_t icase) {
    return m_excit_cases[icase];
}

void ExcitGenGroup::set_probs(const v_t<prob_t>& probs) {
    DEBUG_ASSERT_EQ(probs.size(), ncase(), "incorrect number of probabilities given");
    prob_t norm = std::accumulate(probs.cbegin(), probs.cend(), 0.0);
    DEBUG_ASSERT_GE(norm, 1e-8, "prob vector norm is too small");
    m_probs = probs;
    for (auto& prob : m_probs) {
        DEBUG_ASSERT_GE(prob, 0.0, "no element of the (un-)normalized prob vector can be negative");
        prob/=norm;
    }
    update_cumprobs();
}

prob_t ExcitGenGroup::get_prob(uint_t icase) const {
    DEBUG_ASSERT_LT(icase, ncase(), "excit gen case index OOB");
    return m_probs[icase];
}

const v_t<prob_t>& ExcitGenGroup::get_probs() const {
    return m_probs;
}

void ExcitGenGroup::log() const {
    v_t<strv_t> rows = {{"Excitation Signature", "Description", "Probability"}};
    for (uint_t icase=0ul; icase<ncase(); ++icase){
        auto exsig_str = exsig::to_string(m_excit_cases[icase].m_exsig);
        auto prob_str = convert::to_string(get_prob(icase));
        rows.push_back({exsig_str, m_excit_cases[icase].m_excit_gen->m_description, prob_str});
    }
    log::info_table("Excitation generation breakdown", rows, true);
}

void ExcitGenGroup::set_probs(const sys::Particles& particles) {
    v_t<prob_t> probs;
    probs.reserve(ncase());
    for (const auto& excase: m_excit_cases)
        probs.push_back(excase.m_excit_gen->approx_nconn(excase.m_exsig, particles));
    set_probs(probs);
}
