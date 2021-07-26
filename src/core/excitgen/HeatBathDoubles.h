//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include <src/core/excitgen/FermionExcitationGenerator.h>
#include "src/core/sample/Aliaser.h"
#include "src/core/field/Fields.h"

/*
 * precomputed sampler for doubles
 */

class HeatBathDoubles : public FermionExcitationGenerator {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian *h, PRNG &prng);

    bool draw(const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    bool draw(const fields::FrmBosOnv &src_onv, fields::FrmBosOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
        return draw(src_onv.m_frm, dst_onv.m_frm, occs, vacs, prob, helem, conn.m_frm);
    }

};

#endif //M7_HEATBATHDOUBLES_H
