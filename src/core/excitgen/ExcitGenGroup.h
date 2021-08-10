//
// Created by rja on 04/08/2021.
//

#ifndef M7_EXCITGENGROUP_H
#define M7_EXCITGENGROUP_H

#include "src/core/parallel/Reduction.h"
#include "Hubbard1dSingles.h"
#include "HeatBathDoubles.h"

/**
 * dynamically constructs an array of those excitation generators required for the stochastic propagation of the given
 * Hamiltonian. All excitation signatures which in general give rise to non-zero H matrix elements are called "active"
 * exsigs, and are conventionally indexed with the symbol iex.
 * E.g.
 * The ab-initio fermionic many-body Hamiltonian has two such active exsigs: 0 (1100), 1 (2200)
 * whereas the 1d site-basis Hubbard-Holstein model has three active exsigs: 0 (1100), 1 (0010), 2 (0001)
 */
class ExcitGenGroup {
    PRNG &m_prng;
    /**
     * one excitation generator is dynamically constructable per active exsig
     */
    std::array<std::unique_ptr<ExcitGen>, defs::nexsig> m_exgens;
    /**
     * vector storing all active exsigs consecutively
     */
    defs::inds m_active_exsigs;
    /**
     * probability of attempting to draw from each of the active exsigs
     */
    std::vector<defs::prob_t> m_probs;
    /**
     * cached cumulative probability for all the active exsigs for slight performance benefit
     */
    std::vector<defs::prob_t> m_cumprobs;

    /**
     * indices which point to positions in m_active_exsigs corresponding to purely fermionic excitations
     */
    defs::inds m_frm_inds;
    /**
     * probabilities for each purely fermionic active exsig
     */
    std::vector<defs::prob_t> m_frm_probs;
    std::vector<defs::prob_t> m_frm_cumprobs;
    /**
     * normalization for the fermionic active exsigs
     */
    defs::prob_t m_frm_norm = 1.0;

    /**
     * initialize vectors of exsigs and pointers to excitation generators in general. Also initialize the vector
     * m_frm_inds, which identifies positions of purely fermionic excitation generators from the general vector
     */
    void init();

    /**
     * given only the current value of m_probs, update the probability of purely fermionic exctiation generations and
     * update the cumulative probability vectors for all excitations and for fermionic excitations only
     */
    void update_cumprobs();

public:
    /**
     * infer from the given Hamiltonian exactly which excitation generators are required, and initialize them
     * @param ham
     *  general Hamiltonian object whose excitation level information is queried to determine the required exlvl-specific
     *  excitation generation objects
     * @param prng
     *  random number generator needed to construct the required exlvl-specific excitation generators and decide which
     *  exlvl to attempt to draw
     */
    ExcitGenGroup(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng);

    /**
     * @result
     *  the number of active exsigs
     */
    size_t size() const;

    void set_probs(const std::vector<defs::prob_t> &probs);

    defs::prob_t get_prob(const size_t &iex) const;

    defs::prob_t get_prob_frm(const size_t &iex) const;

    const std::vector<defs::prob_t> &get_probs() const;

    ExcitGen &operator[](const size_t &iex);

    const ExcitGen &operator[](const size_t &iex) const;

    size_t draw_iex();

    size_t draw_iex_frm();

    void log_breakdown() const;
};


#endif //M7_EXCITGENGROUP_H
