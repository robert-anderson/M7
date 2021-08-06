//
// Created by rja on 04/08/2021.
//

#ifndef M7_EXCITGENGROUP_H
#define M7_EXCITGENGROUP_H

#include "src/core/parallel/Reduction.h"
#include "Hubbard1dSingles.h"
#include "HeatBathDoubles.h"


class ExcitGenGroup {
    PRNG &m_prng;
    std::unique_ptr<FrmExcitGen> m_frm_singles;
    std::unique_ptr<FrmExcitGen> m_frm_doubles;
    std::unique_ptr<FrmBosExcitGen> m_frmbos;
    std::vector<defs::prob_t> m_probs;
    std::vector<defs::prob_t> m_cumprobs;
    std::vector<ExcitGen *> m_ptrs;

    void init_probs();

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

    size_t size() const;

    void set_probs(const std::vector<defs::prob_t> &probs);

private:
    void set_prob(const defs::prob_t &prob, ExcitGen *ptr);

public:

    const defs::prob_t &get_prob(const size_t &iexlvl) const;

    const std::vector<defs::prob_t> &get_probs() const;

    ExcitGen &operator[](const size_t &iexlvl);

    size_t draw_iexlvl();

    void log_breakdown() const;
};


#endif //M7_EXCITGENGROUP_H
