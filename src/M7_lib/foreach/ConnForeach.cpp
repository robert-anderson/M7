//
// Created by rja on 28/03/2022.
//

#include "ConnForeach.h"

conn_foreach::Base::Base(size_t exsig, const BasisData &bd) : m_exsig(exsig), m_bd(bd), m_conns(bd){}


void conn_foreach::Base::loop(conn::FrmOnv &conn, const field::FrmOnv& src, const conn_foreach::Base::function_t<conn::FrmOnv> &fn) {
    frm_loop(conn, src, fn);
}

void conn_foreach::Base::loop(const field::FrmOnv& src, const conn_foreach::Base::function_t<conn::FrmOnv> &fn) {
    loop(m_conns.m_frmonv, src, fn);
}

void conn_foreach::Base::loop(conn::BosOnv &conn, const field::BosOnv& src, const conn_foreach::Base::function_t<conn::BosOnv> &fn) {
    bos_loop(conn, src, fn);
}

void conn_foreach::Base::loop(const field::BosOnv& src, const conn_foreach::Base::function_t<conn::BosOnv> &fn) {
    loop(m_conns.m_bosonv, src, fn);
}

void conn_foreach::Base::loop(conn::FrmBosOnv &conn, const field::FrmBosOnv& src, const conn_foreach::Base::function_t<conn::FrmBosOnv> &fn) {
    frmbos_loop(conn, src, fn);
}

void conn_foreach::Base::loop(const field::FrmBosOnv& src, const conn_foreach::Base::function_t<conn::FrmBosOnv> &fn) {
    loop(m_conns.m_frmbosonv, src, fn);
}

void conn_foreach::frm::Base::frmbos_loop(conn::FrmBosOnv &conn, const field::FrmBosOnv &src,
                                          const function_t<conn::FrmBosOnv> &fn) {
    auto frm_fn = [&conn, &fn](const conn::FrmOnv& frm_conn) {
        fn(conn);
    };
    frm_loop(conn.m_frm, src.m_frm, frm_fn);
}
