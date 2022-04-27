//
// Created by Robert J. Anderson on 03/04/2022.
//

#include "ExcitGen.h"

bool ExcitGen::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, conn::FrmBosOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::FrmOnv &conn) {
    auto success = draw_h_frm(exsig, src, prob, helem, conn);
    return success &! consts::nearly_zero(helem, defs::helem_tol);
}

bool ExcitGen::draw(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::FrmBosOnv &conn) {
    auto success = draw_h_frmbos(exsig, src, prob, helem, conn);
    return success &! consts::nearly_zero(helem, defs::helem_tol);
}

bool ExcitGen::draw(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                    conn::BosOnv &conn) {
    auto success = draw_h_bos(exsig, src, prob, helem, conn);
    return success &! consts::nearly_zero(helem, defs::helem_tol);
}
