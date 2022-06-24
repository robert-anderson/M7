//
// Created by Robert J. Anderson on 06/04/2022.
//

#ifndef M7_EXCITGENGROUP_H
#define M7_EXCITGENGROUP_H

#include "M7_lib/hamiltonian/Hamiltonian.h"

using namespace exsig;

/**
 * excitation generator objects can (if uncommonly) be used to generate excitations of many different signatures. e.g.
 * a fermion-boson excitation generator from which an excitation of exsigs 1101 or 1110 may be requested. This struct
 * packages the polymorphic excitation generator base class pointer with one of the exsigs it is able to produce
 */
struct ExcitCase {
    const uint_t m_exsig;
    ExcitGen* m_excit_gen;
};

class ExcitGenGroup {
    /**
     * random number generator used to select which excitation case to attempt
     */
    PRNG& m_prng;
    /**
     * forward list containing all excitation generators for all hamiltonian terms
     */
    ExcitGen::excit_gen_list_t m_list;
    /**
     * one dynamically-constructed excitation generator can be used to generate excitations of many different exsigs
     */
    std::vector<ExcitCase> m_excit_cases;
    /**
     * mapping from the excitation signature to a vector case indices that can generate it
     */
    std::vector<uintv_t> m_exsig_icases;
    /**
     * probability of attempting to draw from each of the active excitation cases
     */
    std::vector<prob_t> m_probs;
    /**
     * cached cumulative probability for all the active excitation cases for slight performance benefit
     */
    std::vector<prob_t> m_cumprobs;

    void update_cumprobs();

public:
    ExcitGenGroup(const Hamiltonian &ham, const conf::Propagator &opts, PRNG &prng);
    /*
     * ctor which also sets heterogeneous initial probabilities to cases by calling the approx_nconn method of exgens
     */
    ExcitGenGroup(const Hamiltonian &ham, const conf::Propagator &opts, PRNG &prng, sys::Particles particles):
            ExcitGenGroup(ham, opts, prng) {
        set_probs(particles);
    }

    uint_t ncase() const;

    uint_t draw_icase();
    
    ExcitCase& operator[](uint_t icase);

    /**
     * set the m_probs and m_cumprobs arrays
     * @param probs
     */
    void set_probs(const std::vector<prob_t> &probs);

    /**
     * set probs in accordance with the approximate number of connections of each excitation case.
     * useful for initialization
     * @param particles
     *  particle number data on which to base the approximate number of connections
     */
    void set_probs(const sys::Particles& particles);

    prob_t get_prob(uint_t icase) const;

    const std::vector<prob_t>& get_probs() const;

    /**
     * when there is strictly one excitation generator per exsig, the probability of drawing the connection is
     * p(conn|exsig) p(exsig)
     * in general, however there's more than one way to draw the same connection
     * p(conn) = sum_case p(conn|case) p(case)
     * this is the update performed here.
     * Naturally if there is only one case for the given exsig, the probability is scaled by the probability of the case.
     */
    template<typename mbf_t>
    void update_prob(uint_t icase, const mbf_t &src, prob_t &prob, const conn::from_field_t<mbf_t> &conn) {
        const auto exsig = m_excit_cases[icase].m_exsig;
        DEBUG_ASSERT_EQ(exsig, conn.exsig(), "exsig of case does not match with that of connection");
        const auto& jcases = m_exsig_icases[exsig];
        DEBUG_ASSERT_FALSE(jcases.empty(), "there should be at least one case associated with this exsig");
        prob *= m_probs[icase];
        if (jcases.size()==1) return;
        for (const auto& jcase: jcases) {
            if (jcase==icase) continue; // already included above
            prob += m_excit_cases[icase].m_excit_gen->prob(src, conn)*m_probs[icase];
        }
    }

    template<typename mbf_t>
    bool draw(uint_t icase, const mbf_t &src, prob_t &prob, ham_t &helem, conn::from_field_t<mbf_t> &conn) {
        auto& excase = m_excit_cases[icase];
        return excase.m_excit_gen->draw(excase.m_exsig, src, prob, helem, conn);
    }

    void log() const;

};


#endif //M7_EXCITGENGROUP_H
