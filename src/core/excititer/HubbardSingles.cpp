//
// Created by rja on 24/11/2021.
//

#include "HubbardSingles.h"

excititers::HubbardSingles::HubbardSingles(const Hamiltonian &ham) :
        Frm(ham, exsig_utils::ex_single), m_pbc(true) {
//    REQUIRE_TRUE(ham.m_frm.is_hubbard_1d() || ham.m_frm.is_hubbard_1d_pbc(),
//                 "Hamiltonian is not 1D hubbard, so this class will not generate all connections");
}

void excititers::HubbardSingles::foreach(const FrmOnv &src, conn::FrmOnv &conn,
                                           const fn_c_t<FrmOnv> &body) {
    const auto &occs = m_work_orbs.occ(src).m_flat.inds();
    for (const auto &occ: occs) {
        size_t neighbor;
        neighbor = model_utils::left(occ, m_bd.m_nsite, m_pbc);
        if (neighbor != ~0ul && !src.get(neighbor)) {
            conn.set(occ, neighbor);
            set_helement(src, conn);
            body(conn);
        }
        neighbor = model_utils::right(occ, m_bd.m_nsite, m_pbc);
        if (neighbor != ~0ul && !src.get(neighbor)) {
            conn.set(occ, neighbor);
            set_helement(src, conn);
            body(conn);
        }
    }
}