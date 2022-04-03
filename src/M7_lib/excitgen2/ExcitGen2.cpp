//
// Created by rja on 03/04/2022.
//

#include "ExcitGen2.h"

bool ExcitGen2::draw_frm(const size_t &exsig, const field::FrmOnv &src, defs::prob_t &prob, conn::FrmOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen2::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src, defs::prob_t &prob, conn::FrmBosOnv &conn) {
    prob = 0.0;
    return false;
}

bool ExcitGen2::draw_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, conn::BosOnv &conn) {
    prob = 0.0;
    return false;
}