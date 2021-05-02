//
// Created by rja on 02/05/2021.
//

#ifndef M7_HUBBARDSINGLES_H
#define M7_HUBBARDSINGLES_H

#include "UniformSingles.h"

struct HubbardSingles : public UniformSingles {
    const bool m_pbc;

    HubbardSingles(const Hamiltonian<>* h, PRNG& prng, bool pbc):
        UniformSingles(h, prng), m_pbc(pbc){}

    bool _draw(const fields::Onv<0> &src_onv, fields::Onv<0> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<0> &anticonn);

    bool _draw(const fields::Onv<1> &src_onv, fields::Onv<1> &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::Antisym<1> &anticonn) {
        return _draw(src_onv.m_frm, dst_onv.m_frm, occs, vacs, prob, helem, anticonn);
    }

public:
    bool draw(const fields::Onv<> &src_onv, fields::Onv<> &dst_onv, const OccupiedOrbitals &occs,
              const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
              conn::Antisym<> &anticonn) override {
        return _draw(src_onv, dst_onv, occs, vacs, prob, helem, anticonn);
    }

};


#endif //M7_HUBBARDSINGLES_H
