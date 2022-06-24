//
// Created by rja on 24/11/2021.
//

#include "LadderPure.h"


excititers::LadderPure::LadderPure(const Hamiltonian &ham, uint_t exsig) : Ladder(ham, exsig) {}

void excititers::LadderPure::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    conn.clear();
    for (uint_t imode=0ul; imode<src.m_bos.m_nelement; ++imode){
        if (m_cre) {
            if (src.m_bos[imode] == m_ham.m_nboson_max) continue;
            conn.m_bos.m_cre.set(imode);
        }
        else {
            if (src.m_bos[imode] == 0ul) continue;
            conn.m_bos.m_ann.set(imode);
        }
        if (!set_helement(src, conn)) continue;
        set_helement(src, conn);
        body(conn);
    }
}