//
// Created by rja on 29/02/2020.
//

#ifndef M7_DETERMINANTSAMPLER_H
#define M7_DETERMINANTSAMPLER_H


#include "HeatBathSampler.h"
#include "HeatBathExcitation.h"


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
    std::vector<defs::prob_t> m_P3;
    std::vector<defs::prob_t> m_P4;

public:
    DeterminantSampler(const HeatBathSampler &precomputed, const Determinant &det);

    static std::vector<defs::prob_t> make_P1(const HeatBathSampler &precomputed,
                                             const defs::inds &occinds, const size_t &noccind);

    void make_P2(std::vector<defs::prob_t> &P2, const size_t &p);
    void make_P3(std::vector<defs::prob_t> &P3, const size_t &p, const size_t &q);
    void make_P4(std::vector<defs::prob_t> &P4, const size_t &p, const size_t &q, const size_t &r);

    HeatBathExcitation draw(PRNG &prng);

private:
    defs::prob_t proposal(const size_t &ip, const size_t &r, const defs::ham_t &h);

    defs::prob_t proposal(const size_t &ip, const size_t &iq,
                          const size_t &r, const size_t &s,
                          const defs::ham_t &h);
};


#endif //M7_DETERMINANTSAMPLER_H
