//
// Created by Robert J. Anderson on 03/04/2022.
//

#include "FrmExcitGen.h"

#include <utility>

FrmExcitGen::FrmExcitGen(const FrmHam &h, PRNG &prng, uintv_t exsigs, str_t description) :
        ExcitGen(prng, std::move(exsigs), std::move(description)), m_h(h){
    for (auto exsig: m_exsigs)
        REQUIRE_TRUE(exsig::is_pure_frm(exsig), "excitations must be expressed in terms of fermion operators only");
}

bool FrmExcitGen::draw_frmbos(uint_t exsig, const field::FrmBosOnv &src,
                              prob_t &prob, conn::FrmBosOnv &conn) {
    return draw_frm(exsig, src.m_frm, prob, conn.m_frm);
}

bool FrmExcitGen::draw_h_frm(uint_t exsig, const field::FrmOnv &src, prob_t &prob,
                             ham_t &helem, conn::FrmOnv &conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return ham::is_significant(helem);
}

bool FrmExcitGen::draw_h_frmbos(uint_t exsig, const field::FrmBosOnv &src, prob_t &prob,
                                ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src.m_frm, prob, conn.m_frm);
    if (!result) return false;
    helem = m_h.get_element(src.m_frm, conn.m_frm);
    return ham::is_significant(helem);
}