//
// Created by rja on 24/11/2021.
//

#include "HubbardSingles.h"

excititers::HubbardSingles::HubbardSingles(const Hamiltonian &ham) :
        Frm(ham, exsig_utils::ex_single), m_pbc(true) {
}

void excititers::HubbardSingles::foreach(const FrmOnv &src, conn::FrmOnv &conn,
                                           const fn_c_t<FrmOnv> &body) {
    const auto &occs = m_work_orbs.occ(src).m_flat.inds();
    auto h_cast = dynamic_cast<const HubbardFrmHam*>(m_ham.m_frm.get());
    REQUIRE_TRUE(h_cast, "Fermion hamiltonian is not hubbard type");
    auto& sparse_conns = h_cast->m_t_mat_sparse;
    for (auto &occ: occs) {
        for (auto& neighbor: sparse_conns[occ].first) {
            conn.set(occ, neighbor);
            set_helement(src, conn);
            body(conn);
        }
    }
}