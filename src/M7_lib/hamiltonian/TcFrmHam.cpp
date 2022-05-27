//
// Created by anderson on 01/04/2022.
//

#include "TcFrmHam.h"

defs::ham_t TcFrmHam::get_coeff_3300(size_t a, size_t b, size_t c, size_t i,
                                     size_t j, size_t k) const {
    auto ia = FrmOnvField::isite(a, m_nsite);
    auto ib = FrmOnvField::isite(b, m_nsite);
    auto ic = FrmOnvField::isite(c, m_nsite);
    auto ii = FrmOnvField::isite(i, m_nsite);
    auto ij = FrmOnvField::isite(j, m_nsite);
    auto ik = FrmOnvField::isite(k, m_nsite);
    // check spin conservation
    auto ma = FrmOnvField::ispin(a, m_nsite);
    auto mb = FrmOnvField::ispin(b, m_nsite);
    auto mc = FrmOnvField::ispin(c, m_nsite);
    auto mi = FrmOnvField::ispin(i, m_nsite);
    auto mj = FrmOnvField::ispin(j, m_nsite);
    auto mk = FrmOnvField::ispin(k, m_nsite);
    defs::ham_t coeff = 0.0;
    auto add_contrib = [&](size_t mp, size_t mq, size_t mr, size_t ip,
                           size_t iq, size_t ir, int sgn) {
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

defs::ham_t TcFrmHam::get_element_3300(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto element = get_coeff_3300(conn.m_cre[0], conn.m_cre[1], conn.m_cre[2],
                                  conn.m_ann[0], conn.m_ann[1], conn.m_ann[2]);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t TcFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    auto element = GeneralFrmHam::get_element_0000(onv);
    auto triples_fn = [&](size_t i, size_t j, size_t k) {
        element += get_coeff_3300(i, j, k, i, j, k);
    };
    onv.foreach_setbit_triple(triples_fn);
    return element;
}

defs::ham_t TcFrmHam::get_element_1100(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto general_element = GeneralFrmHam::get_element_1100(onv, conn);
    defs::ham_t tc_element = 0.0;
    auto doubles_fn = [&](size_t i, size_t j) {
        if (i == conn.m_ann[0] || j == conn.m_ann[0]) return;
        tc_element += get_coeff_3300(conn.m_cre[0], i, j, conn.m_ann[0], i, j);
    };
    onv.foreach_setbit_pair(doubles_fn);
    // inelegant way to get around using the phase twice
    // TODO modify to calculate phase only once
    return general_element + (conn.phase(onv) ? -tc_element : tc_element);
}

defs::ham_t TcFrmHam::get_element_2200(const field::FrmOnv &onv,
                                       const conn::FrmOnv &conn) const {
    auto general_element = GeneralFrmHam::get_element_2200(onv, conn);
    defs::ham_t tc_element = 0.0;
    auto fn = [&](size_t i) {
        if (i == conn.m_ann[0] || i == conn.m_ann[1]) return;
        tc_element += get_coeff_3300(conn.m_cre[0], conn.m_cre[1], i,
                                  conn.m_ann[0], conn.m_ann[1], i);
    };
    onv.foreach_setbit(fn);
    // inelegant way to get around using the phase twice
    // TODO modify to calculate phase only once
    return general_element + (conn.phase(onv) ? -tc_element : tc_element);
}
