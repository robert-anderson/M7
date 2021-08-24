//
// Created by rja on 01/06/2021.
//

#include <src/core/util/Foreach.h>
#include "ForeachConnection.h"
#include "Hamiltonian.h"

foreach_conn::Base::Base(const Hamiltonian &ham) :
        m_ham(ham), m_occ(ham.nsite()), m_vac(ham.nsite()), m_conns(ham.nsite()), m_mbfs(ham.nsite()) {}


defs::ham_t foreach_conn::Base::get_element(const FrmOnv &mbf, const conn::FrmOnv &conn) {
    return m_ham.get_element(mbf, conn);
}

defs::ham_t foreach_conn::Base::get_element(const FrmBosOnv &mbf, const conn::FrmBosOnv &conn) {
    return m_ham.get_element(mbf, conn);
}

defs::ham_t foreach_conn::Base::get_element(const BosOnv &mbf, const conn::BosOnv &conn) {
    return m_ham.get_element(mbf, conn);
}

void foreach_conn::frm::Fermion::singles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
        for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
            conn.set(m_occ[iocc], m_vac[ivac]);
            fn();
        }
    }
}

void foreach_conn::frm::Fermion::singles(
        const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    singles(conn, fn);
}

void foreach_conn::frm::Fermion::doubles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
        for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
            for (size_t jocc = 0ul; jocc < iocc; ++jocc) {
                for (size_t jvac = 0ul; jvac < ivac; ++jvac) {
                    conn.set(m_occ[jocc], m_occ[iocc], m_vac[jvac], m_vac[ivac]);
                    fn();
                }
            }
        }
    }
}

void
foreach_conn::frm::Fermion::foreach(const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    m_occ.update(mbf);
    m_vac.update(mbf);
    singles(conn, fn);
    doubles(conn, fn);
}

void foreach_conn::frm::Fermion::foreach(const FrmOnv &mbf, const fn_c_t<FrmOnv> &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmonv); };
    this->foreach(mbf, m_conns.m_frmonv, fn);
}

void foreach_conn::frm::Fermion::foreach(const FrmBosOnv &mbf, const fn_c_t<FrmBosOnv> &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmbosonv); };
    this->foreach(mbf.m_frm, m_conns.m_frmbosonv.m_frm, fn);
}

void foreach_conn::frm::SpinSym::singles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    Fermion::singles(conn, fn);
}

void foreach_conn::frm::SpinSym::doubles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    Fermion::doubles(conn, fn);
}

foreach_conn::frm::Hubbard1D::Hubbard1D(const Hamiltonian &ham) :
        Fermion(ham), m_pbc(ham.m_frm.m_int_1(0, ham.nsite() - 1) != 0.0) {}

void foreach_conn::frm::Hubbard1D::singles(
        const FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (const auto &occ: m_occ.inds()) {
        size_t neighbor;
        neighbor = model_utils::left(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !mbf.get(neighbor)) {
            conn.set(occ, neighbor);
            fn();
        }
        neighbor = model_utils::right(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !mbf.get(neighbor)) {
            conn.set(occ, neighbor);
            fn();
        }
    }
}

void foreach_conn::frm::Hubbard1D::foreach(const FrmOnv &mbf, conn::FrmOnv &conn,
                                           const std::function<void()> &fn) {
    m_occ.update(mbf);
    singles(mbf, conn, fn);
}



void foreach_conn::bos::Boson::foreach(const FrmBosOnv &mbf, const foreach_conn::fn_c_t<FrmBosOnv> &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmbosonv); };
    this->foreach(mbf.m_bos, m_conns.m_frmbosonv.m_bos, fn);
}

void foreach_conn::bos::Boson::foreach(const BosOnv &mbf, const foreach_conn::fn_c_t<BosOnv> &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_bosonv); };
    this->foreach(mbf, m_conns.m_bosonv, fn);
}

void foreach_conn::bos::Boson::foreach(const BosOnv &mbf, conn::BosOnv &conn, const std::function<void()> &fn) {
    for (size_t imode=0ul; imode<mbf.m_nelement; ++imode){
        if (mbf[imode] < m_ham.m_nboson_max) {
            conn.clear();
            conn.m_cre.set({imode, 1ul});
            DEBUG_ASSERT_EQ(conn.m_bos.size(), 1ul, "should have only one boson operator index");
            fn();
        }
        if (mbf[imode] > 0) {
            conn.clear();
            conn.m_ann.set({imode, 1ul});
            DEBUG_ASSERT_EQ(conn.m_bos.size(), 1ul, "should have only one boson operator index");
            fn();
        }
    }
}


void foreach_conn::frmbos::Holstein::foreach(
        const FrmBosOnv &mbf, const foreach_conn::fn_c_t<FrmBosOnv> &body_fn) {
    if (!m_ham.m_nboson_max) return;
    auto &conn = m_conns.m_frmbosonv;
    for (size_t imode = 0ul; imode < mbf.m_bos.nelement(); ++imode) {
        if (mbf.m_frm.site_nocc(imode)) {
            if (mbf.m_bos[imode] < m_ham.m_nboson_max) {
                conn.clear();
                conn.m_bos.m_cre.set({imode, 1ul});
                DEBUG_ASSERT_EQ(conn.m_bos.size(), 1ul, "should have only one boson operator index");
                body_fn(conn);
            }
            if (mbf.m_bos[imode] > 0) {
                conn.clear();
                conn.m_bos.m_ann.set({imode, 1ul});
                DEBUG_ASSERT_EQ(conn.m_bos.size(), 1ul, "should have only one boson operator index");
                body_fn(conn);
            }
        }
    }
}

void foreach_conn::frmbos::FrmBos::foreach(const FrmBosOnv &mbf, const foreach_conn::fn_c_t<FrmBosOnv> &body_fn) {
    // do all the 0001 / 0010 excits first
    auto &conn = m_boson.m_conns.m_frmbosonv;
    conn.clear();

    m_boson.foreach(mbf, body_fn);

    // now do all those that involve single-fermion excitations
    m_occ.update(mbf);
    m_vac.update(mbf);
    for (const auto& occ : m_occ.inds()){
        for (const auto& vac : m_vac.inds()){
            conn.clear();
            conn.m_frm.set(occ, vac);
            if (!conn.m_frm.kramers_conserve()) continue;
            m_boson.foreach(mbf, body_fn);
        }
    }
}