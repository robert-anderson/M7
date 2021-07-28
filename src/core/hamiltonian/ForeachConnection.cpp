//
// Created by rja on 01/06/2021.
//

#include <src/core/util/Foreach.h>
#include "ForeachConnection.h"
#include "Hamiltonian.h"

foreach_conn::Base::Base(const Hamiltonian &ham) :
        m_ham(ham), m_occ(ham.nsite()), m_vac(ham.nsite()), m_conns(ham.nsite()), m_mbfs(ham.nsite()) {}


defs::ham_t foreach_conn::Base::get_element(const fields::FrmOnv &mbf, const conn::FrmOnv &conn) {
    return m_ham.get_element(mbf, conn);
}

defs::ham_t foreach_conn::Base::get_element(const fields::FrmBosOnv &mbf, const conn::FrmBosOnv &conn) {
    return m_ham.get_element(mbf, conn);
}

defs::ham_t foreach_conn::Base::get_element(const fields::BosOnv &mbf, const conn::BosOnv &conn) {
    return m_ham.get_element(mbf, conn);
}



void foreach_conn::frm::Fermion::singles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
        for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
            conn.clear();
            conn.add(m_occ[iocc], m_vac[ivac]);
            fn();
        }
    }
}

void foreach_conn::frm::Fermion::singles(
        const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    singles(conn, fn);
}

void foreach_conn::frm::Fermion::doubles(conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (size_t iocc = 0ul; iocc < m_occ.size(); ++iocc) {
        for (size_t ivac = 0ul; ivac < m_vac.size(); ++ivac) {
            for (size_t jocc = 0ul; jocc < iocc; ++jocc) {
                for (size_t jvac = 0ul; jvac < ivac; ++jvac) {
                    conn.clear();
                    conn.add(m_occ[jocc], m_occ[iocc], m_vac[jvac], m_vac[ivac]);
                    fn();
                }
            }
        }
    }
}

void
foreach_conn::frm::Fermion::foreach(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    m_occ.update(mbf);
    m_vac.update(mbf);
    singles(conn, fn);
    doubles(conn, fn);
}

void foreach_conn::frm::Fermion::foreach(const fields::FrmOnv &mbf, const foreach_conn::fn_t<fields::FrmOnv> &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmonv); };
    this->foreach(mbf, m_conns.m_frmonv, fn);
}

void foreach_conn::frm::Fermion::foreach(const fields::FrmBosOnv &mbf, const foreach_conn::fn_t<fields::FrmBosOnv> &body_fn) {
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
        const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    for (const auto &occ: m_occ.inds()) {
        size_t neighbor;
        neighbor = conn_utils::left(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !mbf.get(neighbor)) {
            conn.set(occ, neighbor);
            fn();
        }
        neighbor = conn_utils::right(occ, m_ham.nsite(), m_pbc);
        if (neighbor != ~0ul && !mbf.get(neighbor)) {
            conn.set(occ, neighbor);
            fn();
        }
    }
}

void
foreach_conn::frm::Hubbard1D::foreach(const fields::FrmOnv &mbf, conn::FrmOnv &conn,
                                         const std::function<void()> &fn) {
    m_occ.update(mbf);
    singles(mbf, conn, fn);
}
