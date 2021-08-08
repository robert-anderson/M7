//
// Created by rja on 02/05/2021.
//

#ifndef M7_HUBBARD1DSINGLES_H
#define M7_HUBBARD1DSINGLES_H

#include "UniformSingles.h"

struct Hubbard1dSingles : public UniformSingles {
    const bool m_pbc;

    Hubbard1dSingles(const Hamiltonian& h, PRNG& prng, bool pbc):
        UniformSingles(h, prng), m_pbc(pbc){}

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


#endif //M7_HUBBARD1DSINGLES_H
