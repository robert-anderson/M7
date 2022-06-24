//
// Created by Robert J. Anderson on 03/04/2022.
//

#ifndef M7_FRMEXCITGEN_H
#define M7_FRMEXCITGEN_H

#include "M7_lib/excitgen/ExcitGen.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"

struct FrmExcitGen : ExcitGen {

    const FrmHam& m_h;

    FrmExcitGen(const FrmHam& h, PRNG& prng, uintv_t exsigs, std::string description);

    bool draw_frmbos(uint_t exsig, const field::FrmBosOnv& src,
                     prob_t& prob, conn::FrmBosOnv& conn) override;

    bool draw_h_frm(uint_t exsig, const field::FrmOnv& src, prob_t& prob,
                    ham_t& helem, conn::FrmOnv& conn) override;

    bool draw_h_frmbos(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob,
                       ham_t& helem, conn::FrmBosOnv& conn) override;

    bool draw_h_bos(uint_t /*exsig*/, const field::BosOnv& /*src*/, prob_t& prob,
                    ham_t& helem, conn::BosOnv& /*conn*/) override {
        prob = 0.0;
        helem = 0.0;
        return false;
    }

};

struct FrmLatticeExcitGen : FrmExcitGen {
protected:
    mutable lattice::adj_row_t m_work_adj_row;
public:
    FrmLatticeExcitGen(const FrmHam& h, PRNG& prng, uintv_t exsigs, std::string description):
            FrmExcitGen(h, prng, exsigs, description){
        REQUIRE_TRUE(h.m_basis.m_lattice.get(), "Lattice excitation generator requires lattice definition in basis");
    }
};


#endif //M7_FRMEXCITGEN_H
