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

    bool draw(const fields::FrmOnv &src_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) override;

    bool draw(const fields::FrmBosOnv &src_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) override {
        return draw(src_onv.m_frm, occs, vacs, prob, helem, conn.m_frm);
    }

    size_t approx_nconn() const override;

};

#endif //M7_HEATBATHDOUBLES_H
