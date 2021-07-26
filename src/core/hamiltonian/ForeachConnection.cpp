//
// Created by rja on 01/06/2021.
//

#include <src/core/util/Foreach.h>
#include "ForeachConnection.h"
#include "Hamiltonian.h"

foreach_conn::Base::Base(const Hamiltonian &ham) :
        m_ham(ham), m_occ(ham.nsite()), m_vac(ham.nsite()), m_conns(ham.nsite()), m_mbfs(ham.nsite()) {}

void foreach_conn::Base::operator()(const fields::FrmOnv &mbf,
                                    const foreach_conn::frm_h_fn_t &body_fn, bool nonzero_h_only) {
    frm_fn_t fn = [&](const conn::FrmOnv &conn) {
        auto helem = m_ham.get_element(mbf, conn);
        if (nonzero_h_only && consts::float_is_zero(helem)) return;
        body_fn(conn, helem);
    };
    (*this)(mbf, fn);
}

void foreach_conn::Base::operator()(const fields::BosOnv &mbf,
                                    const foreach_conn::bos_h_fn_t &body_fn, bool nonzero_h_only) {
    bos_fn_t fn = [&](const conn::BosOnv &conn) {
        auto helem = m_ham.get_element(mbf, conn);
        if (nonzero_h_only && consts::float_is_zero(helem)) return;
        body_fn(conn, helem);
    };
    (*this)(mbf, fn);
}

void foreach_conn::Base::operator()(const fields::FrmBosOnv &mbf,
                                    const foreach_conn::frmbos_h_fn_t &body_fn, bool nonzero_h_only) {
    frmbos_fn_t fn = [&](const conn::FrmBosOnv &conn) {
        auto helem = m_ham.get_element(mbf, conn);
        if (nonzero_h_only && consts::float_is_zero(helem)) return;
        body_fn(conn, helem);
    };
    (*this)(mbf, fn);
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
foreach_conn::frm::Fermion::operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn, const std::function<void()> &fn) {
    m_occ.update(mbf);
    m_vac.update(mbf);
    singles(conn, fn);
    doubles(conn, fn);
}

void foreach_conn::frm::Fermion::operator()(const fields::FrmOnv &mbf, const foreach_conn::frm_fn_t &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmonv); };
    (*this)(mbf, m_conns.m_frmonv, fn);
}

void foreach_conn::frm::Fermion::operator()(const fields::FrmBosOnv &mbf, const foreach_conn::frmbos_fn_t &body_fn) {
    auto fn = [&]() { body_fn(m_conns.m_frmbosonv); };
    (*this)(mbf.m_frm, m_conns.m_frmbosonv.m_frm, fn);
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
foreach_conn::frm::Hubbard1D::operator()(const fields::FrmOnv &mbf, conn::FrmOnv &conn,
                                         const std::function<void()> &fn) {
    m_occ.update(mbf);
    singles(mbf, conn, fn);
}
