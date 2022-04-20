//
// Created by rja on 05/04/2022.
//

#ifndef M7_PCHB1101HC_H
#define M7_PCHB1101HC_H

#include "FrmBosExcitGen.h"
#include "M7_lib/sample/Aliaser.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"

/**
 * pre-compute a one-per-node shared Aliaser akin to that of the Pchb2200 but this time for sampling the 1101 exsig and
 * its hermitian conjugate 1110 of the general fermion-boson Hamiltonian
 */
struct Pchb1101hc : public FrmBosOpenExcitGen {
    Aliaser m_pick_n_given_pq;

    Pchb1101hc(const FrmBosHam &h, size_t nboson_max, PRNG &prng);

    bool draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, conn::FrmBosOnv &conn) override;
};



#endif //M7_PCHB1101HC_H
