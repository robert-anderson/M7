//
// Created by rja on 04/08/2021.
//

#ifndef M7_EXCITGENGROUP_H
#define M7_EXCITGENGROUP_H

#include "Hubbard1dSingles.h"
#include "HeatBathDoubles.h"


class ExcitGenGroup {
    PRNG& m_prng;
    std::unique_ptr<FrmExcitGen> m_frm_singles;
    std::unique_ptr<FrmExcitGen> m_frm_doubles;
    std::unique_ptr<FrmBosExcitGen> m_frmbos;
    std::vector<defs::prob_t> m_probs;
    std::vector<defs::prob_t> m_cumprobs;
    std::vector<ExcitGen *> m_ptrs;

    defs::inds m_ndraws;

    void init_probs() {
        for (auto ptr: m_ptrs) m_probs.push_back(ptr->approx_nconn());
        auto norm = std::accumulate(m_probs.cbegin(), m_probs.cend(), 0.0);
        for (auto &prob: m_probs) prob /= norm;
        update_cumprobs();
    }

    void update_cumprobs() {
        m_cumprobs = m_probs;
        for(size_t iprob = 1ul; iprob<m_probs.size(); ++iprob)
            m_cumprobs[iprob] = m_cumprobs[iprob-1]+m_probs[iprob];
        DEBUG_ASSERT_TRUE(consts::floats_nearly_equal(1.0, m_cumprobs.back()),
                          "cumulative probability should be 1.0");
    }

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
    ExcitGenGroup(const Hamiltonian &ham, const fciqmc_config::Propagator &opts, PRNG &prng): m_prng(prng) {
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
        m_ndraws.assign(size(), 0ul);
    }

    size_t size() const {
        return m_ptrs.size();
    }

    void set_probs(const std::vector<defs::prob_t>& probs) {
        DEBUG_ASSERT_EQ(probs.size(), size(), "incorrect number of probabilities given");
        m_probs = probs;
        update_cumprobs();
    }

private:
    void set_prob(const defs::prob_t& prob, ExcitGen* ptr){
        DEBUG_ASSERT_LE(prob, 1.0, "given prob exceeds one");
        DEBUG_ASSERT_GE(prob, 0.0, "given prob is less than zero");
        auto norm_remain = 1.0-prob;
        for (size_t i=0ul; i<size(); ++i) {
            if (m_ptrs[i]!=ptr) m_probs[i]*=norm_remain;
        }
        update_cumprobs();
    }

public:
    void set_prob_frm_singles(const defs::prob_t& prob) {
        set_prob(prob, m_frm_singles.get());
    }
    void set_prob_frm_doubles(const defs::prob_t& prob) {
        set_prob(prob, m_frm_doubles.get());
    }
    void set_prob_frmbos(const defs::prob_t& prob) {
        set_prob(prob, m_frmbos.get());
    }

    const defs::prob_t& get_prob(const size_t &i) const {
        DEBUG_ASSERT_LT(i, size(), "excit gen index OOB");
        return m_probs[i];
    }

    ExcitGen &operator[](const size_t &i) {
        DEBUG_ASSERT_LT(i, size(), "excit gen index OOB");
        return *m_ptrs[i];
    }

    size_t draw_exlvl_ind(){
        auto r = m_prng.draw_float();
        for (size_t i=0ul; i < size(); ++i) {
            if (r<m_cumprobs[i]) return i;
        }
        return size()-1;
    }

    bool draw(const fields::FrmOnv &src_onv,
                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn){
        auto exlvl_ind = draw_exlvl_ind();
        auto res = m_ptrs[exlvl_ind]->draw(src_onv, occs, vacs, prob, helem, conn);
        prob*=m_probs[exlvl_ind];
        return res;
    }

    bool draw(const fields::FrmBosOnv &src_onv,
                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn){
        auto exlvl_ind = draw_exlvl_ind();
        auto res = m_ptrs[exlvl_ind]->draw(src_onv, occs, vacs, prob, helem, conn);
        prob*=m_probs[exlvl_ind];
        return res;
    }

    void log_breakdown() const {
        log::info("Excitation class probability breakdown:");
        for (size_t i=0ul; i<size(); ++i)
            log::info("{:<40} {}", m_ptrs[i]->description(), get_prob(i));
    }
};


#endif //M7_EXCITGENGROUP_H
