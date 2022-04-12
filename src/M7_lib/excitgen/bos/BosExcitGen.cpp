//
// Created by rja on 05/04/2022.
//

#include "BosExcitGen.h"

BosExcitGen::BosExcitGen(const BosHam &h, PRNG &prng, defs::inds exsigs, std::string description) :
        ExcitGen(prng, std::move(exsigs), std::move(description)), m_h(h) {
    for (auto exsig: m_exsigs)
        REQUIRE_TRUE(exsig_utils::is_pure_bos(exsig), "excitations must be expressed in terms of boson operators only");
}

bool BosExcitGen::draw_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                              defs::prob_t &prob, conn::FrmBosOnv &conn) {
    return draw_bos(exsig, src.m_bos, prob, conn.m_bos);
}

bool BosExcitGen::draw_h_frm(const size_t &exsig, const field::FrmOnv &src,
                             defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) {
    helem = 0.0;
    return false;
}

bool BosExcitGen::draw_h_frmbos(const size_t &exsig, const field::FrmBosOnv &src,
                                defs::prob_t &prob, defs::ham_t &helem, conn::FrmBosOnv &conn) {
    auto result = draw(exsig, src.m_bos, prob, conn.m_bos);
    if (!result) return false;
    helem = m_h.get_element(src.m_bos, conn.m_bos);
    return !consts::nearly_zero(helem, defs::helem_tol);
}

bool BosExcitGen::draw_h_bos(const size_t &exsig, const field::BosOnv &src, defs::prob_t &prob, defs::ham_t &helem,
                             conn::BosOnv &conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return !consts::nearly_zero(helem, defs::helem_tol);
}
