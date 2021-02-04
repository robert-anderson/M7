//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FermionExcitationGenerator.h"

class UniformSingles : public FermionExcitationGenerator {

public:
    UniformSingles(const Hamiltonian<>* ham, PRNG& prng);

private:
    bool _draw(const views::Onv<0> &src_onv, views::Onv<0> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<0> &anticonn);

    bool _draw(const views::Onv<1> &src_onv, views::Onv<1> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<1> &anticonn) {
        return _draw(src_onv.m_fonv, dst_onv.m_fonv, occs, vacs, prob, helem, anticonn);
    }

public:
    bool draw(const views::Onv<> &src_onv, views::Onv<> &dst_onv, const OccupiedOrbitals &occs,
                    const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
                    conn::Antisym<> &anticonn) override {
        return _draw(src_onv, dst_onv, occs, vacs, prob, helem, anticonn);
    }

};


#endif //M7_UNIFORMSINGLES_H
