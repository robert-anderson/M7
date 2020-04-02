//
// Created by rja on 29/02/2020.
//

#ifndef M7_DETERMINANTSAMPLER_H
#define M7_DETERMINANTSAMPLER_H


#include "HeatBathSampler.h"
#include "HeatBathExcitation.h"

#if 0
class DeterminantSampler {
    const HeatBathSampler &m_precomputed;
    const Determinant &m_det;
    const defs::inds m_occinds;
    const size_t m_noccind;
    const defs::inds m_uncinds;
    const size_t m_nuncind;
    std::vector<defs::prob_t> m_P1;
    const Aliaser m_P1_aliaser;
    std::vector<defs::prob_t> m_P2_qp; // prob of picking q given p picked first
    std::vector<defs::prob_t> m_P2_pq; // prob of picking p given q picked first

public:
    DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det);

    static std::vector<defs::prob_t> make_P1(const HeatBathSampler &precomputed,
                                             const defs::inds &occinds, const size_t &noccind);

    void make_P2(std::vector<defs::prob_t> &P2, const size_t &p);

    void draw_pq(PRNG &prng, size_t &p, size_t &q);

    void draw_pq(PRNG &prng, size_t &p, size_t &q, defs::prob_t &prob);

    void draw_r(PRNG &prng, const size_t &p, const size_t &q, size_t &r);

    void draw_r(PRNG &prng, const size_t &p, const size_t &q, size_t &r, defs::prob_t &prob);

    void draw_pqr(PRNG &prng, size_t &p, size_t &q, size_t &r);

    void draw_pqr(PRNG &prng, size_t &p, size_t &q, size_t &r, defs::prob_t &prob);

    void draw_s(PRNG &prng, const size_t &p, const size_t &q, const size_t &r, size_t &s);

    void draw_s(PRNG &prng, const size_t &p, const size_t &q, const size_t &r, size_t &s, defs::prob_t &prob);

    void draw(PRNG &prng, size_t &p, size_t &q, size_t &r, size_t &s,
              defs::prob_t &prob_single, defs::prob_t &prob_double,
              defs::ham_t &helement_single, defs::ham_t &helement_double);

    HeatBathExcitation draw(PRNG &prng);



private:
    defs::prob_t proposal(const size_t &ip, const size_t &r, const defs::ham_t &helement_single);

    defs::prob_t proposal(const size_t &ip, const size_t &iq,
                          const size_t &r, const size_t &s,
                          const defs::ham_t &h);
};


#endif //M7_DETERMINANTSAMPLER_H
#endif //M7_DETERMINANTSAMPLER_H
