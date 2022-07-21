//
// Created by Robert J. Anderson on 12/9/21.
//

#ifndef M7_HOLSTEINLADDERHAM_H
#define M7_HOLSTEINLADDERHAM_H

#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"

struct HolsteinLadderHam : FrmBosHam {

    const ham_t m_g;

    HolsteinLadderHam(const sys::frm::Basis& frm_basis, ham_t g, uint_t bos_occ_cutoff=sys::bos::c_max_occ) :
        FrmBosHam({frm_basis, {frm_basis.m_nsite, bos_occ_cutoff}}), m_g(g) {
        m_contribs_1101.set_nonzero(exsig::ex_0001);
        m_contribs_1110.set_nonzero(exsig::ex_0010);
    }

    ham_t get_coeff_1110(uint_t imode, uint_t a, uint_t i) const override;

    ham_t get_coeff_1101(uint_t imode, uint_t a, uint_t i) const override;

    ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override;

    excit_gen_list_t make_excit_gens(PRNG& prng, const conf::Propagator& propagator) const override;

    conn_foreach::base_list_t make_foreach_iters() const override;

};


#endif //M7_HOLSTEINLADDERHAM_H
