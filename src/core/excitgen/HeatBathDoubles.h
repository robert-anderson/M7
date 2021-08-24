//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include "ExcitGen.h"
#include "src/core/sample/Aliaser.h"
#include "src/core/field/Fields.h"

/**
 * precomputed sampler for fermion double excitations
 */
class HeatBathDoubles : public FrmExcitGen {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return draw(exsig, src.m_frm, orbs, prob, helem, conn.m_frm);
    }

    size_t approx_nconn() const override;

};

#endif //M7_HEATBATHDOUBLES_H
