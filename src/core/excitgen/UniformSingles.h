//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "FermionExcitationGenerator.h"

class UniformSingles : public FermionExcitationGenerator {

public:
    UniformSingles(const Hamiltonian* ham, PRNG& prng);

    bool draw(const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn);

    bool draw(const fields::FrmBosOnv &src_onv, fields::FrmBosOnv &dst_onv,
               const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
               defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
        return draw(src_onv.m_frm, dst_onv.m_frm, occs, vacs, prob, helem, conn.m_frm);
    }

};


#endif //M7_UNIFORMSINGLES_H
