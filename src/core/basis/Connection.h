//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_CONNECTION_H
#define M7_CONNECTION_H

#include <src/core/field/DeterminantSpecifier.h>
#include <algorithm>
#include <src/core/field/Views.h>

/*
 * <bra| ...cre... ...ann... |ket>
 * the m_ann array contains the spin orbital indices in the ket but not the bra
 * the m_cre array contains the spin orbital indices in the bra but not the ket
 *
 * connect: given two determinants (ket, bra), compute the ann and cre arrays
 * which transforms ket into bra.
 *
 * excite: given a ket determinant, copy its value to a given bra determinant,
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

class Connection {
    const size_t m_element_dsize;
protected:
    const size_t m_nbit;
    defs::det_work m_ann{}, m_cre{};
    size_t m_nann, m_ncre;

public:
    explicit Connection(const DeterminantSpecifier& field);
    Connection(const views::Determinant &ket, const views::Determinant &bra);
    explicit Connection(const views::Determinant &ket);

    const defs::det_work& ann() const {return m_ann;}
    const size_t& ann(const size_t& i) const {return m_ann[i];}
    const size_t& nann() const {return m_nann;}

    const defs::det_work& cre() const {return m_cre;}
    const size_t& cre(const size_t& i) const {return m_cre[i];}
    const size_t& ncre() const {return m_ncre;}

    virtual void connect(const views::Determinant &ket, const views::Determinant &bra);
    virtual void apply(const views::Determinant &ket, views::Determinant &bra);
    void zero(){m_ncre=0; m_nann=0;}
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
class AntisymConnection : public Connection {
    defs::det_work m_com{};
    size_t m_ncom;
    bool m_phase;

public:
    explicit AntisymConnection(const DeterminantSpecifier& field);
    AntisymConnection(const views::Determinant &ket, const views::Determinant &bra);
    explicit AntisymConnection(const views::Determinant &ket);

    void connect(const views::Determinant &ket, const views::Determinant &bra) override;
    void apply(const views::Determinant &ket);
    void apply(const views::Determinant &ket, views::Determinant &bra) override;

    const defs::det_work& com() const {return m_com;}
    //void com(const defs::inds &v) {m_com.assign(v.begin(), v.end()); m_ncom=v.size();}
    const size_t& com(const size_t& i) const {return m_com[i];}
    const size_t& ncom() const {return m_ncom;}
    const bool& phase() const {return m_phase;}
};

template<typename T>
struct MatrixElement {
    AntisymConnection aconn;
    T element = 0;
    MatrixElement(const views::Determinant& det): aconn(AntisymConnection(det)) {}
};

#endif //M7_CONNECTION_H