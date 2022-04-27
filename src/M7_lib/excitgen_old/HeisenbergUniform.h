//
// Created by Robert J. Anderson on 2/6/22.
//

#ifndef M7_HEISENBERGUNIFORM_H
#define M7_HEISENBERGUNIFORM_H

#include "ExcitGen.h"

struct HeisenbergUniform : FrmExcitGen {

    const HeisenbergFrmHam* h_cast() {
        return dynamic_cast<const HeisenbergFrmHam*>(m_h.m_frm.get());
    }
    const HeisenbergFrmHam* h_cast() const {
        return dynamic_cast<const HeisenbergFrmHam*>(m_h.m_frm.get());
    }

    HeisenbergUniform(const Hamiltonian& h, PRNG& prng): FrmExcitGen(h, prng, exsig_utils::ex_double){
        REQUIRE_TRUE(h_cast(), "fermion hamiltonian is not of heisenberg hamiltonian type");
    }

    bool draw_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                  conn::FrmOnv &conn) override;

    std::string description() const override;

    size_t approx_nconn() const override;
};


#endif //M7_HEISENBERGUNIFORM_H
