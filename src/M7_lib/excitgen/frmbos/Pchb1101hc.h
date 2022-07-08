//
// Created by Robert J. Anderson on 05/04/2022.
//

#ifndef M7_PCHB1101HC_H
#define M7_PCHB1101HC_H

#include "FrmBosExcitGen.h"
#include "M7_lib/sample/Aliaser.h"
#include "M7_lib/excitgen/frm/UniformSingles.h"

namespace exgen {
/**
 * pre-compute a one-per-node shared Aliaser akin to that of the Pchb2200 but this time for sampling the 1101 exsig and
 * its hermitian conjugate 1110 of the general fermion-boson Hamiltonian
 */
    struct Pchb1101hc : public FrmBosExcitGen {
        Aliaser m_pick_n_given_pq;

        Pchb1101hc(const FrmBosHam &h, PRNG &prng);

        bool draw_frmbos(uint_t exsig, const field::FrmBosOnv &src, prob_t &prob, conn::FrmBosOnv &conn) override;

        prob_t prob_h_frmbos(const field::FrmBosOnv &src, const conn::FrmBosOnv &conn, ham_t helem) const override;

        prob_t prob_frmbos(const field::FrmBosOnv &src, const conn::FrmBosOnv &conn) const override;

        uint_t approx_nconn(uint_t /*exsig*/, sys::Particles /*particles*/) const override;
    };
}

#endif //M7_PCHB1101HC_H
