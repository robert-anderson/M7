//
// Created by Robert J. Anderson on 06/04/2022.
//

#ifndef M7_EXCITGENGROUP_H
#define M7_EXCITGENGROUP_H

#include "M7_lib/hamiltonian/Hamiltonian.h"

using namespace exsig_utils;

/**
 * excitation generator objects can (if uncommonly) be used to generate excitations of many different signatures. e.g.
 * a fermion-boson excitation generator from which an excitation of exsigs 1101 or 1110 may be requested. This struct
 * packages the polymorphic excitation generator base class pointer with one of the exsigs it is able to produce
 */
struct ExcitCase {
    const size_t m_exsig;
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
     * probability of attempting to draw from each of the active excitation cases
     */
    std::vector<defs::prob_t> m_probs;
    /**
     * cached cumulative probability for all the active excitation cases for slight performance benefit
     */
    std::vector<defs::prob_t> m_cumprobs;

    void update_cumprobs();

public:
    ExcitGenGroup(const Hamiltonian &ham, const conf::Propagator &opts, PRNG &prng);

    size_t ncase() const;

    size_t draw_icase();
    
    ExcitCase& operator[](size_t icase);

    void set_probs(const std::vector<defs::prob_t> &probs);

    defs::prob_t get_prob(const size_t &icase) const;

    const std::vector<defs::prob_t>& get_probs() const;

    bool draw(const size_t &icase, const FrmOnv &src, prob_t &prob, ham_t &helem, conn::FrmOnv &conn);

    bool draw(const size_t &icase, const FrmBosOnv &src, prob_t &prob, ham_t &helem, conn::FrmBosOnv &conn);

    bool draw(const size_t &icase, const BosOnv &src, prob_t &prob, ham_t &helem, conn::BosOnv &conn);

    void log() const;

};


#endif //M7_EXCITGENGROUP_H
