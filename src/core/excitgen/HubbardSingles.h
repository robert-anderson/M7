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

    bool draw(const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    bool draw(const fields::FrmBosOnv &src_onv, fields::FrmBosOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
        return draw(src_onv.m_frm, dst_onv.m_frm, occs, vacs, prob, helem, conn.m_frm);
    }
};


#endif //M7_HUBBARDSINGLES_H
