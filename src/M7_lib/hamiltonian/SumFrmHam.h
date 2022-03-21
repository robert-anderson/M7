//
// Created by rja on 15/03/2022.
//

#ifndef M7_SUMFRMHAM_H
#define M7_SUMFRMHAM_H

#include "FrmHam.h"

/**
 * Implementation of h1 + J*h2 where h1 and h2 are arbitrary subclasses of FrmHam
 * @tparam ham_t1
 *  type of first summand (which is unscaled)
 * @tparam ham_t2
 *  type of second summand (which is scaled)
 */
template<typename ham_t1, typename ham_t2>
struct SumFrmHam : FrmHam {
    ham_t1 m_h1;
    ham_t2 m_h2;
    /**
     * scalar multiple of m_h2
     */
    defs::ham_t m_weight;

    SumFrmHam(ham_t1&& h1, ham_t2&& h2, defs::ham_t weight):
        m_h1(std::move(h1)), m_h2(std::move(h2)), m_weight(weight){}

    defs::ham_t get_coeff_1100(size_t i, size_t j) const override {
        return m_h1.get_coeff_1100(i, j) + m_weight * m_h2.get_coeff_1100(i, j);
    }

    defs::ham_t get_coeff_2200(size_t i, size_t j, size_t k, size_t l) const override {
        return m_h1.get_coeff_2200(i, j) + m_weight * m_h2.get_coeff_2200(i, j);
    }

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override {
        return m_h1.get_element_0000(onv) + m_weight * m_h2.get_element_0000(onv);
    }

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override {
        return m_h1.get_element_1100(onv, conn) + m_weight * m_h2.get_element_1100(onv, conn);
    }

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override{
        return m_h1.get_element_2200(onv, conn) + m_weight * m_h2.get_element_2200(onv, conn);
    }
};


#endif //M7_SUMFRMHAM_H
