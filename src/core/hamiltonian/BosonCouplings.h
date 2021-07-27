//
// Created by rja on 05/11/2020.
//

#ifndef M7_BOSONCOUPLINGS_H
#define M7_BOSONCOUPLINGS_H

#include "src/core/connection/Connections.h"

struct BosonCouplings {
    const size_t m_nboson_cutoff, m_nmode;
    const defs::ham_t m_v;

    BosonCouplings(size_t nmode, size_t nboson_cutoff, defs::ham_t v) :
            m_nboson_cutoff(nboson_cutoff), m_nmode(nmode), m_v(v) {}

    defs::ham_t v(const size_t &p, const size_t &q, const size_t &n) const {
        return (p == q && p == n) ? m_v : 0.0;
    }

    defs::ham_t get_element_1(const fields::FrmBosOnv& onv, const conn::FrmBosOnv &conn) const {
        if (conn.m_bos.size()!=1) return 0.0;
        bool is_ann = conn.m_bos.m_ann.size();

        const auto imode = (conn.m_bos.m_ann.size() ? conn.m_bos.m_ann[0] : conn.m_bos.m_cre[0]).m_imode;
        const auto com = size_t(onv.m_bos[imode])-(is_ann?1:0);

        switch (conn.m_frm.size()) {
            case 0: {
                // fermion ONVs do not differ, so sum over occupied spin orbitals
                defs::ham_t res = 0;
                auto fn = [&](const size_t& ibit){
                    auto isite = ibit % m_nmode;
                    const auto occ_fac = std::sqrt(com + 1);
                    res += v(isite, isite, imode) * occ_fac;
                };
                onv.m_frm.foreach(fn);
            }
            case 2: {
                auto isite = conn.m_frm.m_cre[0] % m_nmode;
                auto jsite = conn.m_frm.m_ann[0] % m_nmode;
                DEBUG_ASSERT_NE(isite, jsite, "single fermion excitation violates spin conservation");
                // TODO: remove the below ASSERT in the general case of fermion-boson interaction
                DEBUG_ASSERT_EQ(v(isite, jsite, imode), 0.0,
                                "Hubbard--Holstein only couples fermions and bosons on the same site")
                return v(isite, jsite, imode);
            }
            default:
                return 0;
        }
    }

};

#endif //M7_BOSONCOUPLINGS_H