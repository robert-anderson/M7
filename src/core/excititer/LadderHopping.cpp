//
// Created by rja on 24/11/2021.
//

#include "LadderHopping.h"



excititers::LadderHopping::LadderHopping(const Hamiltonian &ham, size_t exsig) :
        Ladder(ham, exsig) {}

void excititers::LadderHopping::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    conn.clear();
    const auto &occs = m_work_orbs.occ(src.m_frm).m_flat.inds();
    const auto &vacs = m_work_orbs.vac(src.m_frm).m_flat.inds();

    for (size_t imode = 0ul; imode < m_bd.m_nmode; ++imode) {
        if (m_cre) {
            if (src.m_bos[imode] == m_ham.m_nboson_max) continue;
            conn.m_bos.m_cre.set({imode, 1});
        } else {
            if (src.m_bos[imode] == 0ul) continue;
            conn.m_bos.m_ann.set({imode, 1});
        }
        for (auto &occ: occs) {
            conn.m_frm.m_ann.clear();
            conn.m_frm.m_ann.add(occ);
            for (auto &vac: vacs) {
                conn.m_frm.m_cre.clear();
                conn.m_frm.m_cre.add(vac);
                if (!set_helement(src, conn)) continue;
                body(conn);
            }
        }
    }
}