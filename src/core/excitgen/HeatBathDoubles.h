//
// Created by rja on 09/05/2020.
//

#ifndef M7_HEATBATHDOUBLES_H
#define M7_HEATBATHDOUBLES_H

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include <src/core/excitgen/FermionExcitationGenerator.h>
#include "src/core/sample/Aliaser.h"
#include "src/core/fieldz/FieldsZ.h"

/*
 * precomputed sampler for doubles
 */

class HeatBathDoubles : public FermionExcitationGenerator {
    Aliaser m_pick_ab_given_ij;

public:
    HeatBathDoubles(const Hamiltonian<> *h, PRNG &prng);

    bool _draw(const fieldsz::Onv<0> &src_onv, fieldsz::Onv<0> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<0> &anticonn);

    bool _draw(const fieldsz::Onv<1> &src_onv, fieldsz::Onv<1> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<1> &anticonn) {
        return _draw(src_onv.m_fonv, dst_onv.m_fonv, occs, vacs, prob, helem, anticonn);
    }

    bool draw(const fieldsz::Onv<> &src_onv, fieldsz::Onv<> &dst_onv,
                     const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                     defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<> &anticonn) override {
        return _draw(src_onv, dst_onv, occs, vacs, prob, helem, anticonn);
    }

};

#endif //M7_HEATBATHDOUBLES_H
