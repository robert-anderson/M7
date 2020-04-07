//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_CONNECTION_H
#define M7_CONNECTION_H

#include <src/core/fermion/Determinant.h>
#include <algorithm>

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
    defs::inds m_ann, m_cre;
    size_t m_nann, m_ncre;

public:
    explicit Connection(const Field* field);
    Connection(const DeterminantElement &ket, const DeterminantElement &bra);
    explicit Connection(const DeterminantElement &ket);

    const defs::inds& ann() const {return m_ann;}
    void ann(const defs::inds &v) {m_ann.assign(v.begin(), v.end()); m_nann=v.size();}
    const size_t& ann(const size_t& i) const {return m_ann[i];}
    const size_t& nann() const {return m_nann;}

    const defs::inds& cre() const {return m_cre;}
    void cre(const defs::inds &v) {m_cre.assign(v.begin(), v.end()); m_ncre=v.size();}
    const size_t& cre(const size_t& i) const {return m_cre[i];}
    const size_t& ncre() const {return m_ncre;}

    virtual void connect(const DeterminantElement &ket, const DeterminantElement &bra);
    virtual void apply(const DeterminantElement &ket, DeterminantElement &bra);
    void zero(){m_ncre=0; m_nann=0;}
    void add_cre(const size_t &i){m_cre[m_ncre++] = i;}
    void add_ann(const size_t &i){m_ann[m_nann++] = i;}
    void add(const size_t &ann, const size_t &cre){
        assert(m_nann+1<m_nbit);
        assert(m_ncre+1<m_nbit);
        m_ann[m_nann++] = ann; m_cre[m_ncre++] = cre;
    }
    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2){
        assert(m_nann+2<m_nbit);
        assert(m_ncre+2<m_nbit);
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
    defs::inds m_com;
    size_t m_ncom;
    bool m_phase;

public:
    explicit AntisymConnection(const Field* field);
    AntisymConnection(const DeterminantElement &ket, const DeterminantElement &bra);
    explicit AntisymConnection(const DeterminantElement &ket);

    void connect(const DeterminantElement &ket, const DeterminantElement &bra) override;
    void apply(const DeterminantElement &ket);
    void apply(const DeterminantElement &ket, DeterminantElement &bra) override;

    const defs::inds& com() const {return m_com;}
    void com(const defs::inds &v) {m_com.assign(v.begin(), v.end()); m_ncom=v.size();}
    const size_t& com(const size_t& i) const {return m_com[i];}
    const size_t& ncom() const {return m_ncom;}
    const bool& phase() const {return m_phase;}
};

#endif //M7_CONNECTION_H
