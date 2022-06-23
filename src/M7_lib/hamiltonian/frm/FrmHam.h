//
// Created by Robert J. Anderson on 27/02/2020.
//

#ifndef M7_FRMHAM_H
#define M7_FRMHAM_H

#include <cstddef>

#include "M7_lib/hamiltonian/HamOpTerm.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/io/Options.h"
#include "M7_lib/conf/Conf.h"
#include "M7_lib/table/BufferedFields.h"

#include "M7_lib/hamiltonian/HamiltonianData.h"

/**
 * All interactions between the fermionic parts of MBFs are described in this class.
 * currently, upto three-body fermion number-conserving Hamiltonians are supported:
 * H = sum_a,i T[a,i] a+i + sum_ab,ij U[a,b,i,j] a+b+ji + sum_abc,ijk L[a,b,c,i,j,k] a+b+c+ijk
 * a, b, c are conventionally the creation indices
 * i, j, k are conventionally the annihilation indices
 * T, U, and L are symbols given to the coefficient arrays which parametrise H. Often in model systems, these arrays
 * can be defined by a single scalar which applies only to certain sets of spin orbitals (e.g. based on adjacency), but
 * in other situations (including ab-initio) there are no such simple rules governing the elements of these coefficient
 * arrays and so the GeneralFrmHam derived class which stores the T and U coefficient arrays explicitly is needed.
 * For a case involving non-zero L, see TcFrmHam
 *
 * overview of "get" methods:
 *
 * 1. get_coeff_xx00
 *  these refer to the retrieval of the elements of the H-parametrising arrays (T, U, etc.) where xx00 is a rank
 *  signature
 *
 * 2. get_element_xx00
 *  these refer to "promotions" of excitations, in which Wick-contractions over occupied spin orbitals are computed.
 *  and the sum is multiplied by the +/- 1 Fermi phase due to antisymmetry of the determinant basis. The sums include
 *  coeffs of rank signature yy00 where y >= x e.g.
 *      - the "single-replacement matrix element" takes contributions from T, U and L if defined (but obviously not the
 *        core energy, which has rank 0).
 *      - the "double-replacement contraction" takes contributions only from U, and L.
 */
struct FrmHam : HamOpTerm {
    /**
     * a convenient pair of references to the relevant Hamiltonian section and the Basis configuration section
     */
    typedef HamOpTerm::OptPair<conf::FrmHam> opt_pair_t;
    /**
     * properties of the single particle basis
     */
    const sys::frm::Basis m_basis;
    /**
     * core energy
     */
    defs::ham_t m_e_core = 0.0;
    /**
     * zero or non-zero status of exsig contributions (0000, 1100) to the term of ranksig 1100 i.e. whether the
     * one-electron Hamiltonian has non-zero coefficients corresponding to diagonal elements, and/or single excitations
     */
    ham::TermContribs m_contribs_1100;
    /**
     * zero or non-zero status of exsig contributions (0000, 1100, 2200) to the term of ranksig 2200 i.e. whether the
     * two-electron Hamiltonian has non-zero coefficients corresponding to diagonal elements, and/or single excitations,
     * and/or double excitations
     */
    ham::TermContribs m_contribs_2200;
    /**
     * Time reversal symmetry by term rank in the many-body Hamiltonian
     */
    ham::KramersAttributes m_kramers_attrs;
    /**
     * the only compile time constant depended upon by the Hamiltonian is ENABLE_COMPLEX. The program should still
     * behave properly if a real-valued Hamiltonian is provided to a binary with defs::ham_t compiled as a complex type.
     * this flag indicates whether complex arithmetic is required in order to correctly express the coefficients of the
     * fermion Hamiltonian.
     * An error is thrown if the data source for the coefficients is complex valued but the code is compiled for real-
     * valued hamiltonian matrix elements
     */
    bool m_complex_valued = false;

    FrmHam(const sys::frm::Basis& basis);

private:

    /**
     * workspace for computing connections
     */
    mutable suite::Conns m_work_conn;

public:

    FrmHam(const FrmHam& other): FrmHam(other.m_basis){}

    FrmHam& operator=(const FrmHam& other){return *this;}

	virtual ~FrmHam(){}

    /**
     * coefficient of the "1-body" term
     * @param a
     *  creation spin orbital index
     * @param i
     *  annihilation spin orbital index
     * @return
     *  one-electron coefficient of H given by the 2D array T
     */
    virtual defs::ham_t get_coeff_1100(size_t a, size_t i) const {return 0;}
    /**
     * coefficient of the "2-body" term
     * @param a
     *  creation spin orbital index
     * @param b
     *  creation spin orbital index (b > a)
     * @param i
     *  annihilation spin orbital index
     * @param j
     *  annihilation spin orbital index (j > i)
     * @return
     *  two-electron coefficient of H given by the 4D array U
     */
    virtual defs::ham_t get_coeff_2200(size_t a, size_t b, size_t i, size_t j) const {return 0;}
    /**
     * coefficient of the "3-body" term
     * @param a
     *  creation spin orbital index
     * @param b
     *  creation spin orbital index (b > a)
     * @param c
     *  creation spin orbital index (c > b)
     * @param i
     *  annihilation spin orbital index
     * @param j
     *  annihilation spin orbital index (j > i)
     * @param k
     *  annihilation spin orbital index (k > j)
     * @return
     *  three-electron coefficient of H given by the 6D array L
     */
    virtual defs::ham_t get_coeff_3300(size_t a, size_t b, size_t c, size_t i, size_t j, size_t k) const {return 0;}
    /**
     * @param onv
     *  fermionic occupation number vector (Slater determinant)
     * @return
     *  diagonal matrix element of H in terms of all ranks of coefficients
     */
    virtual defs::ham_t get_element_0000(const field::FrmOnv &onv) const {return 0;}
    /**
     * @param ket
     *  fermionic occupation number vector (Slater determinant)
     * @param conn
     *  represents the normal-ordered product of fermion operators that produces the bra FrmOnv when acted upon ket
     * @return
     *  single-replacement matrix element of H
     */
    virtual defs::ham_t get_element_1100(const field::FrmOnv &ket, const conn::FrmOnv &conn) const {return 0;}
    /**
     * @param ket
     *  fermionic occupation number vector (Slater determinant)
     * @param conn
     *  represents the normal-ordered product of fermion operators that produces the bra FrmOnv when acted upon ket
     * @return
     *  double-replacement matrix element of H
     */
    virtual defs::ham_t get_element_2200(const field::FrmOnv &ket, const conn::FrmOnv &conn) const {return 0;}
    /**
     * @param ket
     *  fermionic occupation number vector (Slater determinant)
     * @param conn
     *  represents the normal-ordered product of fermion operators that produces the bra FrmOnv when acted upon ket
     * @return
     *  triple-replacement matrix element of H
     */
    virtual defs::ham_t get_element_3300(const field::FrmOnv &ket, const conn::FrmOnv &conn) const {return 0;}

    defs::ham_t get_element(const field::FrmOnv &onv) const;

    defs::ham_comp_t get_energy(const field::FrmOnv &onv) const;

    defs::ham_t get_element(const field::FrmOnv &onv, const conn::FrmOnv &conn) const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    virtual void log_data() const;

    virtual size_t default_nelec() const {
        return m_basis.m_nsite;
    }

    virtual int default_ms2_value() const {
        return sys::frm::c_undefined_ms2;
    }
};

/**
 * fermion hamiltonian which can be defined in a non-zero number of sites, but with no non-zero term coefficients
 */
struct NullFrmHam : FrmHam, NullOpTerm {
    NullFrmHam() : FrmHam(0ul){}
};

/**
 * fermion sites may not be doubly occupied or unoccupied in spin systems
 */
struct SpinModelFrmHam : FrmHam {
protected:
    mutable lattice::adj_row_t m_work_adj_row;
public:
    SpinModelFrmHam(const std::shared_ptr<lattice::Base>& lattice): FrmHam(lattice){}
};

#endif //M7_FRMHAM_H
