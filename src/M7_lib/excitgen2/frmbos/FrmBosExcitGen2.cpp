//
// Created by rja on 05/04/2022.
//

#include "FrmBosExcitGen2.h"

FrmBosExcitGen2::FrmBosExcitGen2(const FrmBosHam &h, PRNG &prng, defs::inds exsigs, std::string description) :
        ExcitGen2(prng, std::move(exsigs), std::move(description)), m_h(h) {}

bool FrmBosExcitGen2::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                                 defs::ham_t &helem, conn::FrmOnv &conn) {
    helem = 0.0;
    return false;
}

bool FrmBosExcitGen2::draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                                    defs::ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::nearly_zero(helem, defs::helem_tol);
}

bool FrmBosExcitGen2::draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                                 defs::ham_t &helem, conn::BosOnv &conn) {
    helem = 0.0;
    return false;
}
