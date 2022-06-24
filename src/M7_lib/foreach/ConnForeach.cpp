//
// Created by Robert J. Anderson on 28/03/2022.
//

#include "ConnForeach.h"

conn_foreach::Base::Base(uint_t exsig) : m_exsig(exsig){}

void conn_foreach::Base::loop(conn::FrmOnv& conn, const field::FrmOnv& src, const function_t& fn) {
    frm_loop(conn, src, fn);
}

void conn_foreach::Base::loop(conn::BosOnv& conn, const field::BosOnv& src, const function_t& fn) {
    bos_loop(conn, src, fn);
}

void conn_foreach::Base::loop(conn::FrmBosOnv& conn, const field::FrmBosOnv& src, const function_t& fn) {
    frmbos_loop(conn, src, fn);
}

void conn_foreach::frm::Base::frmbos_loop(conn::FrmBosOnv& conn, const field::FrmBosOnv& src, const function_t& fn) {
    frm_loop(conn.m_frm, src.m_frm, fn);
}
