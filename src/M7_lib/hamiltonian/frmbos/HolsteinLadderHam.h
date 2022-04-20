//
// Created by anderson on 12/9/21.
//

#ifndef M7_HOLSTEINLADDERHAM_H
#define M7_HOLSTEINLADDERHAM_H

#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"

struct HolsteinLadderHam : FrmBosHam {

    const defs::ham_t m_g;

    HolsteinLadderHam(size_t nsite, const FrmHam& frm, const BosHam& bos, defs::ham_t g) :
        FrmBosHam({nsite, nsite}, frm, bos), m_g(g) {
        m_contribs_0001.set_nonzero(exsig_utils::ex_0001);
        m_contribs_0010.set_nonzero(exsig_utils::ex_0010);
    }

    defs::ham_t get_coeff_0010(size_t imode) const override;

    defs::ham_t get_coeff_0001(size_t imode) const override;

    defs::ham_t get_coeff_1110(size_t imode, size_t i, size_t j) const override;

    defs::ham_t get_coeff_1101(size_t imode, size_t i, size_t j) const override;

    defs::ham_t get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const override;

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

};


#endif //M7_HOLSTEINLADDERHAM_H
