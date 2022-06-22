//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "FrmBosExcitGen.h"

FrmBosExcitGen::FrmBosExcitGen(const FrmBosHam &h, PRNG &prng, defs::inds_t exsigs, std::string description) :
        ExcitGen(prng, std::move(exsigs), std::move(description)), m_h(h) {}

bool FrmBosExcitGen::draw_h_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob,
                                defs::ham_t &helem, conn::FrmOnv &conn) {
    helem = 0.0;
    return false;
}

bool FrmBosExcitGen::draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob,
                                   defs::ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::nearly_zero(helem, defs::helem_tol);
}

bool FrmBosExcitGen::draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob,
                                defs::ham_t &helem, conn::BosOnv &conn) {
    helem = 0.0;
    return false;
}