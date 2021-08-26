//
// Created by rja on 25/08/2021.
//

#include "ExcitIters.h"


fn_c_t<FrmOnv>
excititers::FrmConserve::convert(conn::FrmBosOnv &work_conn, const fn_c_t<FrmBosOnv> &fn) {
    work_conn.m_bos.clear();
    return [&](const conn::FrmOnv &conn){
        fn(work_conn);
    };
}

excititers::FrmConserve::FrmConserve(const Hamiltonian &ham, size_t exsig) :
        Frm(ham, exsig),
        m_cre_loop(2*ham.nsite()-ham.nelec(), exsig_utils::decode_nfrm_cre(exsig)),
        m_ann_loop(ham.nelec(), exsig_utils::decode_nfrm_ann(exsig)){
    DEBUG_ASSERT_TRUE(exsig_utils::is_pure_frm(exsig), "exsig should have no bosonic operators");
}

void excititers::FrmConserve::foreach(const FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<FrmOnv> &body) {
    conn.clear();
    auto inner = [&](){
        conn.m_ann.set(m_work_orbs.occ(src).m_flat.inds(), m_ann_loop.inds());
        if (!set_helement(src, conn)) return;
        body(conn);
    };
    auto outer = [&](){
        conn.m_cre.set(m_work_orbs.vac(src).m_flat.inds(), m_cre_loop.inds());
        m_ann_loop(inner);
    };
    m_cre_loop(outer);
}

void excititers::FrmConserve::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn,
                                      const fn_c_t<FrmBosOnv> &body) {
    auto converted_body = convert(conn, body);
    foreach(src.m_frm, conn.m_frm, converted_body);
}



excititers::Hubbard1dSingles::Hubbard1dSingles(const Hamiltonian &ham) :
        FrmConserve(ham, exsig_utils::ex_single), m_pbc(ham.m_frm.is_hubbard_1d_pbc()) {
    REQUIRE_TRUE(ham.m_frm.is_hubbard_1d(),
                 "Hamiltonian is not 1D hubbard, so this class will not generate all connections");
}

void excititers::Hubbard1dSingles::foreach(const FrmOnv &src, conn::FrmOnv &conn,
                                           const fn_c_t<FrmOnv> &body) {
    const auto &occs = m_work_orbs.occ(src).m_flat.inds();
    for (const auto &occ: occs) {
        size_t neighbor;
        neighbor = model_utils::left(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !src.get(neighbor)) {
            conn.set(occ, neighbor);
            set_helement(src, conn);
            body(conn);
        }
        neighbor = model_utils::right(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !src.get(neighbor)) {
            conn.set(occ, neighbor);
            set_helement(src, conn);
            body(conn);
        }
    }
}