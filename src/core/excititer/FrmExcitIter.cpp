//
// Created by rja on 25/08/2021.
//

#include "FrmExcitIter.h"

ExcitIter::fn_c_t<FrmOnv>
FrmExcitIter::convert(conn::FrmBosOnv &work_conn, const ExcitIter::fn_c_t<FrmBosOnv> &fn) {
    work_conn.m_bos.clear();
    return [&](const conn::FrmOnv &conn){
        fn(work_conn);
    };
}

FrmExcitIter::FrmExcitIter(size_t exsig, const Hamiltonian &ham) :
        ExcitIter(exsig, ham),
        m_cre_loop(2*ham.nsite()-ham.nelec(), exsig_utils::decode_nfrm_cre(exsig)),
        m_ann_loop(ham.nelec(), exsig_utils::decode_nfrm_ann(exsig)){
    DEBUG_ASSERT_TRUE(exsig_utils::is_pure_frm(exsig), "exsig should have no bosonic operators");
}

void FrmExcitIter::foreach(const FrmOnv &src, conn::FrmOnv &work_conn, const ExcitIter::fn_c_t<FrmOnv> &body) {
    work_conn.clear();
    auto inner = [&](){
        work_conn.m_ann.set(m_work_orbs.occ(src).m_flat.inds(), m_ann_loop.inds());
        if (!set_helement(src, work_conn)) return;
        body(work_conn);
    };
    auto outer = [&](){
        work_conn.m_cre.set(m_work_orbs.vac(src).m_flat.inds(), m_cre_loop.inds());
        m_ann_loop(inner);
    };
    m_cre_loop(outer);
}

void FrmExcitIter::foreach(const FrmBosOnv &src, conn::FrmBosOnv &work_conn,
                           const ExcitIter::fn_c_t<FrmBosOnv> &body) {
    auto converted_body = convert(work_conn, body);
    foreach(src.m_frm, work_conn.m_frm, converted_body);
}
