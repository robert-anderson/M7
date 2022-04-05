//
// Created by rja on 25/08/2021.
//

#include "ExcitIters.h"


fn_c_t<FrmOnv>
excititers::Frm::convert(conn::FrmBosOnv &work_conn, const fn_c_t<FrmBosOnv> &fn) {
    work_conn.m_bos.clear();
    return [&](const conn::FrmOnv &conn){
        fn(work_conn);
    };
}

excititers::Frm::Frm(const Hamiltonian &ham, size_t exsig) :
        ExcitIter(ham, exsig),
        m_cre_loop(m_bd.m_nspinorb-ham.nelec(), exsig_utils::decode_nfrm_cre(exsig)),
        m_ann_loop(ham.nelec(), exsig_utils::decode_nfrm_ann(exsig)){
    REQUIRE_TRUE(exsig_utils::is_pure_frm(exsig), "excitation signature should not have bosonic operators")
}

void excititers::Frm::foreach(const FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<FrmOnv> &body) {
    conn.clear();
    auto inner = [&](){
        //conn.m_ann.set(m_work_orbs.occ(src).m_flat.inds(), m_ann_loop.inds());
        if (!set_helement(src, conn)) return;
        body(conn);
    };
    auto outer = [&](){
        //conn.m_cre.set(m_work_orbs.vac(src).m_flat.inds(), m_cre_loop.inds());
        m_ann_loop(inner);
    };
    m_cre_loop(outer);
}

void excititers::Frm::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn,
                              const fn_c_t<FrmBosOnv> &body) {
    auto converted_body = convert(conn, body);
    foreach(src.m_frm, conn.m_frm, converted_body);
}

excititers::Ladder::Ladder(const Hamiltonian &ham, size_t exsig) :
        ExcitIter(ham, exsig), m_cre(exsig_utils::decode_nbos_cre(exsig)) {}


fn_c_t<field::BosOnv> excititers::Bos::convert(conn::FrmBosOnv &work_conn, const fn_c_t<FrmBosOnv> &fn) {
    work_conn.m_frm.clear();
    return [&](const conn::BosOnv &conn){
        fn(work_conn);
    };
}

void excititers::Bos::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    auto converted_body = convert(conn, body);
    foreach(src.m_bos, conn.m_bos, converted_body);
}

void excititers::Bos::foreach(const BosOnv &src, conn::BosOnv &conn, const fn_c_t<BosOnv> &body) {
    conn.clear();
    for (size_t imode = 0ul; imode < src.m_size; ++imode) {
        for (size_t jmode = imode; jmode < src.m_size; ++jmode) {
            conn.m_cre.set(imode, jmode);
            for (size_t kmode = 0ul; kmode < src.m_size; ++kmode) {
                if (!src[kmode]) continue;
                if (kmode==imode || kmode==jmode) continue;
                for (size_t lmode = kmode; lmode < src.m_size; ++lmode) {
                    if (!src[lmode]) continue;
                    if (kmode==lmode && src[lmode]==1) continue;
                    if (lmode==imode || lmode==jmode) continue;
                    conn.m_ann.set(kmode, lmode);
                    if (!set_helement(src, conn)) continue;
                    body(conn);
                }
            }
        }
    }
}
