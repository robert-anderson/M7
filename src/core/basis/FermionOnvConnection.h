//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_FERMIONONVCONNECTION_H
#define M7_FERMIONONVCONNECTION_H

#include <src/core/field/FermionOnvSpecifier.h>
#include <algorithm>
#include <src/core/field/Views.h>

/*
 * <out| ...cre... ...ann... |in>
 * the m_ann array contains the spin orbital indices in "in" but not "out"
 * the m_cre array contains the spin orbital indices in "out" but not "in"
 *
 * connect: given two fermion ONVs ("in", "out"), compute the ann and cre arrays
 * which transforms "in" into "out".
 *
 * excite: given a "in" fermion ONV, copy its value to a given "out" fermion ONV,
 * then apply the internal ann and cre arrays to it.
 *
 * zero: reset the counters
 *
 * add_cre: append creation index
 *
 * add_ann: append annihilation index
 *
 * add: append annihilation/creation index pair
 */

class FermionOnvConnection {
    const size_t m_element_dsize;
protected:
    const size_t m_nbit;
    defs::det_work m_ann{}, m_cre{};
    size_t m_nann, m_ncre;

public:
    explicit FermionOnvConnection(const FermionOnvSpecifier& field);
    FermionOnvConnection(const views::FermionOnv &in, const views::FermionOnv &out);
    explicit FermionOnvConnection(const views::FermionOnv &in);

    const defs::det_work& ann() const {return m_ann;}
    const size_t& ann(const size_t& i) const {return m_ann[i];}
    const size_t& nann() const {return m_nann;}

    const defs::det_work& cre() const {return m_cre;}
    const size_t& cre(const size_t& i) const {return m_cre[i];}
    const size_t& ncre() const {return m_ncre;}

    virtual void connect(const views::FermionOnv &in, const views::FermionOnv &out);
    virtual void apply(const views::FermionOnv &in, views::FermionOnv &out);
    virtual void zero(){m_ncre=0; m_nann=0;}
    void add_cre(const size_t &i){m_cre[m_ncre++] = i;}
    void add_ann(const size_t &i){m_ann[m_nann++] = i;}
    void add(const size_t &ann, const size_t &cre){
        ASSERT(m_nann + 1 < m_nbit);
        ASSERT(m_ncre + 1 < m_nbit);
        m_ann[m_nann++] = ann; m_cre[m_ncre++] = cre;
    }
    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2){
        ASSERT(m_nann + 2 < m_nbit);
        ASSERT(m_ncre + 2 < m_nbit);
        m_ann[m_nann++] = ann1; m_ann[m_nann++] = ann2;
        m_cre[m_ncre++] = cre1; m_cre[m_ncre++] = cre2;
    }

    void sort(){
        std::sort(m_cre.begin(), m_cre.begin()+m_ncre);
        std::sort(m_ann.begin(), m_ann.begin()+m_nann);
    }

    const size_t &nexcit() const;

};

/*
 * a connection in which the common indices and antisymmetric phase is computed
 */
class AntisymFermionOnvConnection : public FermionOnvConnection {
    defs::det_work m_com{};
    size_t m_ncom;
    bool m_phase;

public:
    explicit AntisymFermionOnvConnection(const FermionOnvSpecifier& field);
    AntisymFermionOnvConnection(const views::FermionOnv &in, const views::FermionOnv &out);
    explicit AntisymFermionOnvConnection(const views::FermionOnv &in);

    void connect(const views::FermionOnv &in, const views::FermionOnv &out) override;
    void apply(const views::FermionOnv &in);
    void apply(const views::FermionOnv &in, views::FermionOnv &out) override;

    const defs::det_work& com() const {return m_com;}
    const size_t& com(const size_t& i) const {return m_com[i];}
    const size_t& ncom() const {return m_ncom;}
    const bool& phase() const {return m_phase;}

    void zero() override;
};

template<typename T>
struct MatrixElement {
    AntisymFermionOnvConnection aconn;
    T element = 0;
    MatrixElement(const views::FermionOnv& det): aconn(AntisymFermionOnvConnection(det)) {}
};

#endif //M7_FERMIONONVCONNECTION_H
