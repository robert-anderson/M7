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
struct HeatBathDoubles : public FrmExcitGen {
    using FrmExcitGen::draw;
private:
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian &h, PRNG &prng);

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, conn::FrmOnv &conn) override {
        /*
         * need the helement to compute the probability so if it isn't actually needed, just dispose of it
         * in contrast to the generic case where it is not assumed that the helement must be computed to get the prob
         */
        defs::ham_t helem;
        return draw(exsig, src, orbs, prob, helem, conn);
    }

    size_t approx_nconn() const override;

};

#endif //M7_HEATBATHDOUBLES_H
