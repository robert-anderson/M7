//
// Created by anderson on 01/04/2022.
//

#include "TcFrmHam.h"

ham_t TcFrmHam::get_coeff_3300(uint_t a, uint_t b, uint_t c, uint_t i, uint_t j, uint_t k) const {
    auto ia = m_basis.isite(a);
    auto ib = m_basis.isite(b);
    auto ic = m_basis.isite(c);
    auto ii = m_basis.isite(i);
    auto ij = m_basis.isite(j);
    auto ik = m_basis.isite(k);
    // check spin conservation
    auto ma = m_basis.ispin(a);
    auto mb = m_basis.ispin(b);
    auto mc = m_basis.ispin(c);
    auto mi = m_basis.ispin(i);
    auto mj = m_basis.ispin(j);
    auto mk = m_basis.ispin(k);
    ham_t coeff = 0.0;
    auto add_contrib = [&](uint_t mp, uint_t mq, uint_t mr, uint_t ip,
                           uint_t iq, uint_t ir, int sgn) {
        if (mp == ma && mq == mb && mr == mc) {
            coeff += sgn * get_lmat_coeff(ia, ib, ic, ip, iq, ir);
        }
    };
    // add all exchange terms
    // note this is basically copied from the TCHInt code
    add_contrib(mi, mj, mk, ii, ij, ik, 1);
    add_contrib(mj, mk, mi, ij, ik, ii, 1);
    add_contrib(mk, mi, mj, ik, ii, ij, 1);
    add_contrib(mj, mi, mk, ij, ii, ik, -1);
    add_contrib(mi, mk, mj, ii, ik, ij, -1);
    add_contrib(mk, mj, mi, ik, ij, ii, -1);

    return coeff;
}

ham_t TcFrmHam::get_element_3300(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto element = get_coeff_3300(conn.m_cre[0], conn.m_cre[1], conn.m_cre[2],
                                  conn.m_ann[0], conn.m_ann[1], conn.m_ann[2]);
    return conn.phase(onv) ? -element : element;
}

ham_t TcFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    auto element = GeneralFrmHam::get_element_0000(onv);
    auto triples_fn = [&](uint_t i, uint_t j, uint_t k) {
        element += get_coeff_3300(i, j, k, i, j, k);
    };
    onv.foreach_setbit_triple(triples_fn);
    return element;
}

ham_t TcFrmHam::get_element_1100(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto general_element = GeneralFrmHam::get_element_1100(onv, conn);
    ham_t tc_element = 0.0;
    auto doubles_fn = [&](uint_t i, uint_t j) {
        if (i == conn.m_ann[0] || j == conn.m_ann[0]) return;
        tc_element += get_coeff_3300(conn.m_cre[0], i, j, conn.m_ann[0], i, j);
    };
    onv.foreach_setbit_pair(doubles_fn);
    // inelegant way to get around using the phase twice
    // TODO modify to calculate phase only once
    return general_element + (conn.phase(onv) ? -tc_element : tc_element);
}

ham_t TcFrmHam::get_element_2200(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto general_element = GeneralFrmHam::get_element_2200(onv, conn);
    ham_t tc_element = 0.0;
    auto fn = [&](uint_t i) {
        if (i == conn.m_ann[0] || i == conn.m_ann[1]) return;
        tc_element += get_coeff_3300(conn.m_cre[0], conn.m_cre[1], i,
                                  conn.m_ann[0], conn.m_ann[1], i);
    };
    onv.foreach_setbit(fn);
    // inelegant way to get around using the phase twice
    // TODO modify to calculate phase only once
    return general_element + (conn.phase(onv) ? -tc_element : tc_element);
}
