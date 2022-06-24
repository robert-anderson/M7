//
// Created by Robert J. Anderson on 12/9/21.
//

#ifndef M7_HOLSTEINLADDERHAM_H
#define M7_HOLSTEINLADDERHAM_H

#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

struct HolsteinLadderHam : FrmBosHam {

    const defs::ham_t m_g;

    HolsteinLadderHam(const sys::frm::Basis& frm_basis, defs::ham_t g, uint_t bos_occ_cutoff=sys::bos::c_max_occ) :
        FrmBosHam({frm_basis, {frm_basis.m_nsite, bos_occ_cutoff}}), m_g(g) {
        m_contribs_1101.set_nonzero(exsig::ex_0001);
        m_contribs_1110.set_nonzero(exsig::ex_0010);
    }

    defs::ham_t get_coeff_1110(uint_t imode, uint_t i, uint_t j) const override;

    defs::ham_t get_coeff_1101(uint_t imode, uint_t i, uint_t j) const override;

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

};


#endif //M7_HOLSTEINLADDERHAM_H
