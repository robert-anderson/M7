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
protected:
    defs::inds m_ann, m_cre;

public:
    explicit FermionOnvConnection(const FermionOnvSpecifier& field);
    FermionOnvConnection(const views::Onv<0> &in, const views::Onv<0> &out);
    explicit FermionOnvConnection(const views::Onv<0> &in);
    virtual operator bool() const {
        return nexcit();
    }

    const defs::inds& ann() const {return m_ann;}
    const size_t& ann(const size_t& i) const {return m_ann[i];}
    size_t nann() const {return m_ann.size();}

    const defs::inds& cre() const {return m_cre;}
    const size_t& cre(const size_t& i) const {return m_cre[i];}
    size_t ncre() const {return m_cre.size();}

    void connect(const views::Onv<0> &in, const views::Onv<0> &out);
    void apply(const views::Onv<0> &in, views::Onv<0> &out);
    void zero(){
        m_cre.clear();
        m_ann.clear();
    }

    void add_cre(const size_t &i){
        ASSERT(m_cre.size()<m_cre.capacity());
        m_cre.push_back(i);
    }

    void add_ann(const size_t &i){
        ASSERT(m_ann.size()<m_ann.capacity());
        m_ann.push_back(i);
    }

    void add(const size_t &ann, const size_t &cre){
        add_ann(ann);
        add_cre(cre);
    }

    void add(const size_t &ann1, const size_t &ann2, const size_t &cre1, const size_t &cre2){
        add_ann(ann1);
        add_ann(ann2);
        add_cre(cre1);
        add_cre(cre2);
    }

    void sort(){
        std::sort(m_cre.begin(), m_cre.end());
        std::sort(m_ann.begin(), m_ann.end());
    }

    size_t nexcit() const;

    bool connected() const {
        return nexcit()<=2;
    }

};

/*
 * a connection in which the common indices and antisymmetric phase is computed
 */
class AntisymFermionOnvConnection : public FermionOnvConnection {
    defs::inds m_com;
    bool m_phase;

public:
    explicit AntisymFermionOnvConnection(const FermionOnvSpecifier& spec);
    AntisymFermionOnvConnection(const views::Onv<0> &in, const views::Onv<0> &out);
    explicit AntisymFermionOnvConnection(const views::Onv<0> &in);

    void connect(const views::Onv<0> &in, const views::Onv<0> &out);
    void apply(const views::Onv<0> &in);
    void apply(const views::Onv<0> &in, views::Onv<0> &out);

    void zero();
    const defs::inds & com() const {return m_com;}
    const size_t& com(const size_t& i) const {return m_com[i];}
    size_t ncom() const {return m_com.size();}
    const bool& phase() const {return m_phase;}
};

template<typename T>
struct MatrixElement {
    AntisymFermionOnvConnection aconn;
    T element = 0;
    MatrixElement(const views::Onv<0>& det): aconn(AntisymFermionOnvConnection(det)) {}
};

#endif //M7_FERMIONONVCONNECTION_H
