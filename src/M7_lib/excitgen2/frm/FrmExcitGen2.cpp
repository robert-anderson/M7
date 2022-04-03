//
// Created by rja on 03/04/2022.
//

#include "FrmExcitGen2.h"

FrmExcitGen2::FrmExcitGen2(const FrmHam &h, PRNG &prng, defs::inds exsigs, std::string description) :
        ExcitGen2(prng, std::move(exsigs), std::move(description)),
        m_h(h), m_nelec_pair(integer_utils::nspair(h.m_nelec)),
        m_nspinorb_pair(integer_utils::nspair(2*h.m_nsite)) {}

bool FrmExcitGen2::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                               defs::prob_t &prob, conn::FrmBosOnv &conn) {
    return draw_frm(exsig, src.m_frm, prob, conn.m_frm);
}

bool FrmExcitGen2::draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) {
    return ExcitGen2::draw_bos(exsig, src, prob, conn);
}

bool FrmExcitGen2::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                              defs::ham_t &helem, conn::FrmOnv &conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::nearly_zero(helem);
}

bool FrmExcitGen2::draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                                 defs::ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src.m_frm, prob, conn.m_frm);
    if (!result) return false;
    helem = m_h.get_element(src.m_frm, conn.m_frm);
    return !consts::nearly_zero(helem);
}

bool FrmExcitGen2::draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                              defs::ham_t &helem, conn::BosOnv &conn) {
    helem = 0.0;
    return false;
}
