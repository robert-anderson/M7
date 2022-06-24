//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "BosExcitGen.h"

BosExcitGen::BosExcitGen(const BosHam& h, PRNG& prng, defs::ivec_t exsigs, std::string description) :
        ExcitGen(prng, std::move(exsigs), std::move(description)), m_h(h) {
    for (auto exsig: m_exsigs)
        REQUIRE_TRUE(exsig::is_pure_bos(exsig), "excitations must be expressed in terms of boson operators only");
}

bool BosExcitGen::draw_frmbos(size_t exsig, const field::FrmBosOnv& src,
                              defs::prob_t& prob, conn::FrmBosOnv& conn) {
    return draw_bos(exsig, src.m_bos, prob, conn.m_bos);
}

bool BosExcitGen::draw_h_frmbos(size_t exsig, const field::FrmBosOnv& src,
                                defs::prob_t& prob, defs::ham_t& helem, conn::FrmBosOnv& conn) {
    auto result = draw(exsig, src.m_bos, prob, conn.m_bos);
    if (!result) return false;
    helem = m_h.get_element(src.m_bos, conn.m_bos);
    return ham::is_significant(helem);
}

bool BosExcitGen::draw_h_bos(size_t exsig, const field::BosOnv& src, defs::prob_t& prob, defs::ham_t& helem,
                             conn::BosOnv& conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return ham::is_significant(helem);
}
