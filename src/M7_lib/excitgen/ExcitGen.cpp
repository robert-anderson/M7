//
// Created by Robert J. Anderson on 03/04/2022.
//

#include "ExcitGen.h"
#include "M7_lib/hamiltonian/HamiltonianData.h"

bool ExcitGen::draw_frm(uint_t, const field::FrmOnv&, prob_t& prob, conn::FrmOnv& ) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_frmbos(uint_t, const field::FrmBosOnv&, prob_t& prob, conn::FrmBosOnv& ) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw_bos(uint_t, const field::BosOnv&, prob_t& prob, conn::BosOnv& ) {
    prob = 0.0;
    return false;
}

bool ExcitGen::draw(uint_t exsig, const field::FrmOnv& src, prob_t& prob, ham_t& helem,
                    conn::FrmOnv& conn) {
    auto success = draw_h_frm(exsig, src, prob, helem, conn);
    return success && ham::is_significant(helem);
}

bool ExcitGen::draw(uint_t exsig, const field::FrmBosOnv& src, prob_t& prob, ham_t& helem,
                    conn::FrmBosOnv& conn) {
    auto success = draw_h_frmbos(exsig, src, prob, helem, conn);
    return success && ham::is_significant(helem);
}

bool ExcitGen::draw(uint_t exsig, const field::BosOnv& src, prob_t& prob, ham_t& helem,
                    conn::BosOnv& conn) {
    auto success = draw_h_bos(exsig, src, prob, helem, conn);
    return success && ham::is_significant(helem);
}


prob_t ExcitGen::prob_h_frm(const field::FrmOnv& src, const conn::FrmOnv& conn, ham_t /*helem*/) const  {
    return prob_frm(src, conn);
}

prob_t ExcitGen::prob_h_bos(const field::BosOnv& src, const conn::BosOnv& conn, ham_t /*helem*/) const  {
    return prob_bos(src, conn);
}

prob_t ExcitGen::prob_h_frmbos(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn, ham_t /*helem*/) const  {
    return prob_frmbos(src, conn);
}


prob_t ExcitGen::prob(const field::FrmOnv& src, const conn::FrmOnv& conn) const {
    return prob_frm(src, conn);
}

prob_t ExcitGen::prob(const field::BosOnv& src, const conn::BosOnv& conn)  const {
    return prob_bos(src, conn);
}

prob_t ExcitGen::prob(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn) const  {
    return prob_frmbos(src, conn);
}

prob_t ExcitGen::prob(const field::FrmOnv& src, const conn::FrmOnv& conn, ham_t helem) const  {
    return prob_h_frm(src, conn, helem);
}

prob_t ExcitGen::prob(const field::BosOnv& src, const conn::BosOnv& conn, ham_t helem) const  {
    return prob_h_bos(src, conn, helem);
}

prob_t ExcitGen::prob(const field::FrmBosOnv& src, const conn::FrmBosOnv& conn, ham_t helem) const  {
    return prob_h_frmbos(src, conn, helem);
}
