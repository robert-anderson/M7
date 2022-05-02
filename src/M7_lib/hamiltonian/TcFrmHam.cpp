//
// Created by anderson on 01/04/2022.
//

#include "TcFrmHam.h"

defs::ham_t TcFrmHam::get_coeff_3300(size_t a, size_t b, size_t c, size_t i, size_t j, size_t k) const {
    const int ia=a, ib=b, ic=c, ii=i, ij=j, ik=k;
    // the function below returns a float
    return three_body_coeff(&ia, &ib, &ic, &ii, &ij, &ik);
}

defs::ham_t TcFrmHam::get_element_3300(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    auto element = get_coeff_3300(conn.m_cre[0], conn.m_cre[1], conn.m_cre[2],
                                  conn.m_ann[0], conn.m_ann[1], conn.m_ann[2]);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t TcFrmHam::get_element_0000(const field::FrmOnv &onv) const {
    return GeneralFrmHam::get_element_0000(onv);
    // TODO stub requires a sum over bit triples
}

defs::ham_t TcFrmHam::get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    auto element = GeneralFrmHam::get_element_1100(onv, conn);
    auto doubles_fn = [&](size_t i, size_t j) {
        if (i==conn.m_ann[0] || j==conn.m_ann[0]) return;
        element += get_coeff_3300(conn.m_cre[0], i, j, conn.m_ann[0], i, j);
    };
    onv.foreach_setbit_pair(doubles_fn);
    return conn.phase(onv) ? -element : element;
}

defs::ham_t TcFrmHam::get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
    auto element = GeneralFrmHam::get_element_2200(onv, conn);
    auto fn = [&](size_t i) {
        if (i==conn.m_ann[0] || i==conn.m_ann[1]) return;
        element+=get_coeff_3300(conn.m_cre[0], conn.m_cre[1], i,conn.m_ann[0], conn.m_ann[1], i);
    };
    onv.foreach_setbit(fn);
    return conn.phase(onv) ? -element : element;
}
