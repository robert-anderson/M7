//
// Created by anderson on 12/9/21.
//

#ifndef M7_HOLSTEINLADDERHAM_H
#define M7_HOLSTEINLADDERHAM_H

#include "LadderHam.h"

struct HolsteinLadderHam : LadderHam {

    const defs::ham_t m_g;

    HolsteinLadderHam(size_t nsite, size_t nboson_max, defs::ham_t g):
        LadderHam({nsite, nsite}, nboson_max), m_g(g){
        m_contribs_0001.set_nonzero(exsig_utils::ex_0001);
        m_contribs_0010.set_nonzero(exsig_utils::ex_0010);
    }

    defs::ham_t get_coeff_0010(const size_t &imode) const override {
        return 0;
    }

    defs::ham_t get_coeff_0001(const size_t &imode) const override {
        return 0;
    }

    defs::ham_t get_coeff_1110(const size_t &imode, const size_t &j, const size_t &i) const override {
        return m_g;
    }

    defs::ham_t get_coeff_1101(const size_t &imode, const size_t &j, const size_t &i) const override {
        return m_g;
    }

    defs::ham_t get_element_0010(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }

    defs::ham_t get_element_0001(const field::BosOnv &onv, const conn::BosOnv &conn) const override {
        return 0;
    }

    defs::ham_t get_element_0010(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override {
        auto imode = conn.m_bos.m_cre[0].m_imode;
        if (!onv.m_frm.get({0, imode}) && !onv.m_frm.get({1, imode})) return 0;
        return m_g*conn.m_bos.occ_fac(onv.m_bos);
    }

    defs::ham_t get_element_0001(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override {
        auto imode = conn.m_bos.m_ann[0].m_imode;
        if (!onv.m_frm.get({0, imode}) && !onv.m_frm.get({1, imode})) return 0;
        return m_g*conn.m_bos.occ_fac(onv.m_bos);
    }

    defs::ham_t get_element_1110(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override {
        return 0;
    }

    defs::ham_t get_element_1101(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const override {
        return 0;
    }

};


#endif //M7_HOLSTEINLADDERHAM_H
