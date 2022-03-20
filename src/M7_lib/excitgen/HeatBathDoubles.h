//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include <M7_lib/sample/Aliaser.h>
#include <M7_lib/field/Fields.h>

#include "ExcitGen.h"

/**
 * precomputed sampler for fermion double excitations
 */
struct HeatBathDoubles : public FrmExcitGen {
    using FrmExcitGen::draw;
private:
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian &h, PRNG &prng);

    bool draw_frm(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
                  defs::prob_t &prob, conn::FrmOnv &conn) override {
        /*
         * need the helement to compute the probability so if it isn't actually needed, just dispose of it
         * in contrast to the generic case where it is not assumed that the helement must be computed to get the prob,
         * so the delegation between the two virtual methods is reversed with respect to the generic case.
         */
        defs::ham_t helem;
        return draw_h_frm(exsig, src, orbs, prob, helem, conn);
    }

    bool draw_h_frm(const size_t &exsig, const field::FrmOnv &src, CachedOrbs &orbs,
              defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;


    size_t approx_nconn() const override;

    std::string description() const override;

};

#endif //M7_HEATBATHDOUBLES_H
