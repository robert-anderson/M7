//
// Created by rja on 24/11/2021.
//

#include "LadderPureHolstein.h"


excititers::LadderPureHolstein::LadderPureHolstein(const Hamiltonian &ham, size_t exsig) : Ladder(ham, exsig){}

void excititers::LadderPureHolstein::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    const auto &occs = m_work_orbs.occ(src.m_frm).m_flat.inds();
    for (const auto& occ: occs) {
        conn.clear();
        auto imode = src.m_frm.isite(occ);
        // skip if we reach an alpha orbital whose corresponding beta will produce a call to body
        if (occ < src.nsite() && src.m_frm.get({1, imode})) continue;
        if (m_cre) {
            if (src.m_bos[imode] == m_ham.m_nboson_max) continue;
            conn.m_bos.m_cre.set({imode, 1});
        }
        else {
            if (src.m_bos[imode] == 0ul) continue;
            conn.m_bos.m_ann.set({imode, 1});
        }
        set_helement(src, conn);
        body(conn);
    }
}
