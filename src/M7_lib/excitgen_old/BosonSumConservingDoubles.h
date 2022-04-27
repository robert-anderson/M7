//
// Created by Robert J. Anderson on 25/11/2021.
//

#ifndef M7_BOSONSUMCONSERVINGDOUBLES_H
#define M7_BOSONSUMCONSERVINGDOUBLES_H

#include "ExcitGen.h"

struct BosonSumConservingDoubles : public BosExcitGen {
    const size_t m_nboson_pair;
    BosonSumConservingDoubles(const Hamiltonian &h, PRNG &prng);

    bool draw_frmbos(const size_t &exsig, const FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_bos, orbs, prob, conn.m_bos);
    }

    bool draw_bos(const size_t &exsig, const BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob, conn::BosOnv &conn) override;



    std::string description() const override;

};


#endif //M7_BOSONSUMCONSERVINGDOUBLES_H
