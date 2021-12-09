//
// Created by rja on 04/06/2020.
//

#include "ExcitGen.h"

ExcitGen::ExcitGen(const Hamiltonian &h, PRNG &prng, defs::inds exsigs) :
        m_h(h), m_prng(prng),
        m_bd(m_h.m_bd),
        m_nelec(m_h.nelec()),
        m_norb_pair(integer_utils::nspair(m_bd.m_nspinorb)),
        m_nelec_pair(integer_utils::nspair(m_nelec)), m_exsigs(std::move(exsigs)) {}

bool ExcitGen::draw_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs,
                        defs::prob_t &prob, conn::FrmOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                           conn::FrmBosOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_bos(const size_t &exsig, const BosOnv &src, CachedOrbs &orbs,
                        defs::prob_t &prob, conn::BosOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_h_frm(const size_t &exsig, const FrmOnv &src, CachedOrbs &orbs, defs::prob_t &prob, defs::ham_t &helem,
                          conn::FrmOnv &conn) {
    auto result = draw(exsig, src, orbs, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::float_is_zero(helem);
}

bool ExcitGen::draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                             defs::ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src, orbs, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::float_is_zero(helem);
}

bool ExcitGen::draw_h_bos(const size_t &exsig, const field::BosOnv &src, CachedOrbs &orbs, defs::prob_t &prob,
                          defs::ham_t &helem, conn::BosOnv &conn) {
    auto result = draw(exsig, src, orbs, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::float_is_zero(helem);
}


FrmExcitGen::FrmExcitGen(const Hamiltonian &h, PRNG &prng, size_t exsig) :
        ExcitGen(h, prng, {exsig}),
        m_spin_conserving(exsig == exsig_utils::ex_single ?
                          h.m_frm->m_kramers_attrs.m_conserving_singles : h.m_frm->m_kramers_attrs.m_conserving_double) {
    DEBUG_ASSERT_TRUE(exsig_utils::is_pure_frm(exsig),
                      "fermion excitation generator is not applicable to exsigs with bosonic operators");
    DEBUG_ASSERT_EQ(exsig_utils::decode_nfrm_cre(exsig), exsig_utils::decode_nfrm_ann(exsig),
                    "fermion number non-conservation is not currently supported");
}

BosExcitGen::BosExcitGen(const Hamiltonian &h, PRNG &prng, size_t exsig) :
        ExcitGen(h, prng, {exsig}) {
    DEBUG_ASSERT_TRUE(exsig_utils::is_pure_bos(exsig),
                      "boson excitation generator is not applicable to exsigs with fermionic operators");
}

size_t BosExcitGen::approx_nconn() const {
    return m_h.m_bos->m_nboson * 2;
}
