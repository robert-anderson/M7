//
// Created by Robert J. Anderson on 12/9/21.
//

#include "HolsteinLadderHam.h"
#include "M7_lib/excitgen/frmbos/HolsteinUniform.h"
#include "M7_lib/excitgen/frm/Hubbard.h"

ham_t HolsteinLadderHam::get_coeff_1110(uint_t imode, uint_t a, uint_t i) const {
    if (imode != m_basis.m_frm.isite(a)) return 0;
    if (imode != m_basis.m_frm.isite(i)) return 0;
    return m_g;
}

ham_t HolsteinLadderHam::get_coeff_1101(uint_t imode, uint_t a, uint_t i) const {
    return get_coeff_1110(imode, i, a);
}

ham_t HolsteinLadderHam::get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    const auto imode = conn.m_bos.m_cre[0].m_imode;
    const auto nocc_frm = onv.m_frm.site_nocc(imode);
    const auto occ_fac = onv.m_bos.occ_fac(conn.m_bos);
    return m_g * ham_t(nocc_frm * occ_fac);
}

ham_t HolsteinLadderHam::get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
    const auto imode = conn.m_bos.m_ann[0].m_imode;
    const auto nocc_frm = onv.m_frm.site_nocc(imode);
    const auto occ_fac = onv.m_bos.occ_fac(conn.m_bos);
    return m_g * ham_t(nocc_frm * occ_fac);
}

HamOpTerm::excit_gen_list_t HolsteinLadderHam::make_excit_gens(PRNG& prng, const conf::Propagator& /*propagator*/) const {
    excit_gen_list_t list;
    // boosn creation ("excitation")
    list.emplace_front(new exgen::HolsteinUniform0010(*this, prng));
    // boosn annihilation ("de-excitation")
    list.emplace_front(new exgen::HolsteinUniform0001(*this, prng));
    return list;
}

conn_foreach::base_list_t HolsteinLadderHam::make_foreach_iters() const {
    conn_foreach_list_t list;
    // boosn creation ("excitation")
    list.emplace_front(new conn_foreach::bos::Cre);
    // boosn annihilation ("de-excitation")
    list.emplace_front(new conn_foreach::bos::Ann);
    return list;
}
