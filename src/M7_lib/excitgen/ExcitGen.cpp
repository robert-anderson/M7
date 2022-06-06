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


defs::prob_t ExcitGen::prob_h_frm(const field::FrmOnv &src, const conn::FrmOnv &conn, defs::ham_t helem) const  {
    return prob_frm(src, conn);
}

defs::prob_t ExcitGen::prob_h_bos(const field::BosOnv &src, const conn::BosOnv &conn, defs::ham_t helem) const  {
    return prob_bos(src, conn);
}

defs::prob_t ExcitGen::prob_h_frmbos(const field::FrmBosOnv &src, const conn::FrmBosOnv &conn, defs::ham_t helem) const  {
    return prob_frmbos(src, conn);
}


defs::prob_t ExcitGen::prob(const field::FrmOnv &src, const conn::FrmOnv &conn) const {
    return prob_frm(src, conn);
}

defs::prob_t ExcitGen::prob(const field::BosOnv &src, const conn::BosOnv &conn)  const {
    return prob_bos(src, conn);
}

defs::prob_t ExcitGen::prob(const field::FrmBosOnv &src, const conn::FrmBosOnv &conn) const  {
    return prob_frmbos(src, conn);
}

defs::prob_t ExcitGen::prob(const field::FrmOnv &src, const conn::FrmOnv &conn, defs::ham_t helem) const  {
    return prob_h_frm(src, conn, helem);
}

defs::prob_t ExcitGen::prob(const field::BosOnv &src, const conn::BosOnv &conn, defs::ham_t helem) const  {
    return prob_h_bos(src, conn, helem);
}

defs::prob_t ExcitGen::prob(const field::FrmBosOnv &src, const conn::FrmBosOnv &conn, defs::ham_t helem) const  {
    return prob_h_frmbos(src, conn, helem);
}
