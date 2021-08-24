//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "ExcitGen.h"

class UniformSingles : public FrmExcitGen {

public:
    UniformSingles(const Hamiltonian& ham, PRNG& prng);

    bool draw(const size_t& exsig, const field::FrmOnv &src, CachedOrbs& orbs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const size_t& exsig, const field::FrmBosOnv &src, CachedOrbs& orbs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, helem, conn.m_frm);
    }

    size_t approx_nconn() const override;

};


#endif //M7_UNIFORMSINGLES_H
