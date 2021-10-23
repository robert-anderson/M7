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
        m_cre_loop(m_bd.m_nspinorb-ham.nelec(), exsig_utils::decode_nfrm_cre(exsig)),
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
    REQUIRE_TRUE(ham.m_frm.is_hubbard_1d() || ham.m_frm.is_hubbard_1d_pbc(),
                 "Hamiltonian is not 1D hubbard, so this class will not generate all connections");
}

void excititers::Hubbard1dSingles::foreach(const FrmOnv &src, conn::FrmOnv &conn,
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

excititers::LadderPure::LadderPure(const Hamiltonian &ham, size_t exsig) : Ladder(ham, exsig) {}

void excititers::LadderPure::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    conn.clear();
    for (size_t imode=0ul; imode<src.m_bos.m_nelement; ++imode){
        if (m_cre) {
            if (src.m_bos[imode] == m_ham.m_nboson_max) continue;
            conn.m_bos.m_cre.set({imode, 1});
        }
        else {
            if (src.m_bos[imode] == 0ul) continue;
            conn.m_bos.m_ann.set({imode, 1});
        }
        if (!set_helement(src, conn)) continue;
        set_helement(src, conn);
        body(conn);
    }
}

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

excititers::LadderHopping::LadderHopping(const Hamiltonian &ham, size_t exsig) :
        Ladder(ham, exsig) {}

void excititers::LadderHopping::foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<FrmBosOnv> &body) {
    conn.clear();
    const auto& occs = m_work_orbs.occ(src.m_frm).m_flat.inds();
    const auto& vacs = m_work_orbs.vac(src.m_frm).m_flat.inds();

    for (size_t imode=0ul; imode<m_bd.m_nmode; ++imode) {
        if (m_cre) {
            if (src.m_bos[imode] == m_ham.m_nboson_max) continue;
            conn.m_bos.m_cre.set({imode, 1});
        }
        else {
            if (src.m_bos[imode] == 0ul) continue;
            conn.m_bos.m_ann.set({imode, 1});
        }
        for (auto& occ: occs){
            conn.m_frm.m_ann.clear();
            conn.m_frm.m_ann.add(occ);
            for (auto& vac: vacs) {
                conn.m_frm.m_cre.clear();
                conn.m_frm.m_cre.add(vac);
                if (!set_helement(src, conn)) continue;
                body(conn);
            }
        }
    }
}
