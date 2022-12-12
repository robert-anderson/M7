//
// Created by Robert J. Anderson on 15/03/2022.
//

#ifndef M7_SUMFRMHAM_H
#define M7_SUMFRMHAM_H

#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"

/**
 * Implementation of h1 + J*h2 where h1 and h2 are arbitrary subclasses of FrmHam
 * @tparam ham_t1
 *  type of first summand (which is unscaled)
 * @tparam ham_t2
 *  type of second summand (which is scaled)
 */
template<typename ham_t1, typename ham_t2>
struct SumFrmHam : FrmHam {
    static_assert(std::is_base_of<FrmHam, ham_t1>::value, "first template arg must be derived from FrmHam");
    static_assert(std::is_base_of<FrmHam, ham_t2>::value, "second template arg must be derived from FrmHam");
    ham_t1 m_h1;
    ham_t2 m_h2;
    /**
     * scalar multiple of m_h2
     */
    ham_t m_weight;

private:
    /*
     * convenient base-typed casts
     */
    const FrmHam &m_h1_base;
    const FrmHam &m_h2_base;
public:

    SumFrmHam(ham_t1&& h1, ham_t2&& h2, ham_t weight):
            FrmHam(static_cast<const FrmHam&>(h1).m_basis), m_h1(std::move(h1)), m_h2(std::move(h2)), m_weight(weight),
            m_h1_base(m_h1), m_h2_base(m_h2){
        /*
         * combine attributes of the components
         */
        m_contribs_1100 = ham::TermContribs(m_h1_base.m_contribs_1100, m_h2_base.m_contribs_1100);
        m_contribs_2200 = ham::TermContribs(m_h1_base.m_contribs_2200, m_h2_base.m_contribs_2200);
        m_kramers_attrs = ham::KramersAttributes(m_h1_base.m_kramers_attrs, m_h2_base.m_kramers_attrs);
    }

    ham_t get_coeff_1100(uint_t a, uint_t i) const override {
        return m_h1.get_coeff_1100(a, i) + m_weight * m_h2.get_coeff_1100(a, i);
    }

    ham_t get_coeff_2200(uint_t a, uint_t b, uint_t i, uint_t j) const override {
        return m_h1.get_coeff_2200(a, b, i, j) + m_weight * m_h2.get_coeff_2200(a, b, i, j);
    }

    ham_t get_element_0000(const field::FrmOnv &onv) const override {
        return m_h1.get_element_0000(onv) + m_weight * m_h2.get_element_0000(onv);
    }

    ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override {
        return m_h1.get_element_1100(onv, conn) + m_weight * m_h2.get_element_1100(onv, conn);
    }

    ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override{
        return m_h1.get_element_2200(onv, conn) + m_weight * m_h2.get_element_2200(onv, conn);
    }

private:
    template<typename T>
    static bool any_component_general(const T&) {
        return false;
    }
    /**
     * recurse through components in what is assumed to be a nested SumFrmHam type. If any FrmHam can be dynamically
     * cast to GeneralFrmHam, then the SumFrmHam is deemed to be expressed in general form
     * @tparam ham_u1
     *  type of first summand
     * @tparam ham_u2
     *  type of second summand
     * @param h
     *  sum Hamiltonian
     * @return
     *  true if the sum takes general form
     */
    template<typename ham_u1, typename ham_u2>
    static bool any_component_general(const SumFrmHam<ham_u1, ham_u2>& h) {
        return any_component_general(h.m_h1) || any_component_general(h.m_h2);
    }

    static bool any_component_general(const GeneralFrmHam&) {
        return true;
    }

public:

    excit_gen_list_t make_excit_gens(PRNG &prng, const conf::Propagator &opts) const override {
        if (any_component_general(*this)){
            /*
             * some component of this sum takes the most general form, so let the general excitation generator include
             * all components
             */
            return GeneralFrmHam::make_excit_gens(prng, opts, *this);
        }
        /*
         * otherwise, simply concatenate the two lists
         */
        excit_gen_list_t list;
        list.merge(m_h1_base.make_excit_gens(prng, opts));
        list.merge(m_h2_base.make_excit_gens(prng, opts));
        return list;
    }

    conn_foreach::base_list_t make_foreach_iters() const override {
        conn_foreach::base_list_t list;
        list.merge(m_h1_base.make_foreach_iters());
        list.merge(m_h2_base.make_foreach_iters());
        return list;
    }

    bool is_hermitian() const override {
        return m_h1_base.is_hermitian() && m_h2_base.is_hermitian();
    }

    uint_t default_nelec() const override {
        const auto v1 = m_h1_base.default_nelec();
        const auto v2 = m_h2_base.default_nelec();
        REQUIRE_EQ(v1, v2, "conflicting default electron numbers");
        return v1;
    }

    int default_ms2_value() const override {
        const auto v1 = m_h1_base.default_ms2_value();
        const auto v2 = m_h2_base.default_ms2_value();
        REQUIRE_EQ(v1, v2, "conflicting default Ms2 values");
        return v1;
    }
};


#endif //M7_SUMFRMHAM_H
