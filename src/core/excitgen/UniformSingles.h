//
// Created by RJA on 20/11/2020.
//

#ifndef M7_UNIFORMSINGLES_H
#define M7_UNIFORMSINGLES_H

#include "ExcitGen.h"

class UniformSingles : public FrmExcitGen {

public:
    UniformSingles(const Hamiltonian& ham, PRNG& prng);

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


#endif //M7_UNIFORMSINGLES_H
