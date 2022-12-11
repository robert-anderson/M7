//
// Created by Robert J. Anderson on 05/04/2022.
//

#include "FrmBosExcitGen.h"

FrmBosExcitGen::FrmBosExcitGen(const FrmBosHam& h, PRNG& prng, v_t<OpSig> exsigs, str_t description) :
        ExcitGen(prng, std::move(exsigs), std::move(description)), m_h(h) {}


bool FrmBosExcitGen::draw_h_frmbos(OpSig exsig, const field::FrmBosOnv& src, prob_t& prob,
                                   ham_t& helem, conn::FrmBosOnv& conn) {
    auto result = draw(exsig, src, prob, conn);
    if (!result) return false;
    helem = m_h.get_element(src, conn);
    return ham::is_significant(helem);
}