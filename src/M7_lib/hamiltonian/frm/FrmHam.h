//
// Created by rja on 27/02/2020.
//

#ifndef M7_FRMHAM_H
#define M7_FRMHAM_H

#include <cstddef>

#include "M7_lib/hamiltonian/HamOpTerm.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/io/Options.h"
#include "M7_lib/config/FciqmcConfig.h"
#include "M7_lib/integrals/Integrals_1e.h"
#include "M7_lib/integrals/Integrals_2e.h"
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
 *  these refer to the retrieval of the elements of the H-parametrising arrays (T, U, etc.) where xx00 is a rank signature
 *
 * 2. get_element_xx00
 *  these refer to "promotions" of excitations, in which Wick-contractions over occupied spin orbitals are computed.
 *  and the sum is multiplied by the +/- 1 Fermi phase due to antisymmetry of the determinant basis. The sums include
 *  coeffs of rank signature yy00 where y >= x e.g.
 *      - the "single-replacement matrix element" takes contributions from T, U and L if defined (but obviously not the 
 *      core energy, which has rank 0).
 *      - the "double-replacement contraction" takes contributions only from U, and L.
 */
struct FrmHam : HamOpTerm {

    const size_t m_nelec;
    const size_t m_nsite;
    const int m_ms2_restrict;

    AbelianGroupMap m_point_group_map;

    defs::ham_t m_e_core = 0.0;

    ham_data::TermContribs m_contribs_1100;
    ham_data::TermContribs m_contribs_2200;
    ham_data::KramersAttributes m_kramers_attrs;
    bool m_complex_valued = false;

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict, const defs::inds& site_irreps = {});

    FrmHam(const FrmHam& other): FrmHam(other.m_nelec, other.m_nsite, other.m_ms2_restrict,
                                        other.m_point_group_map.m_site_irreps){}

    FrmHam& operator=(const FrmHam& other){return *this;}

	virtual ~FrmHam(){}

    /**
     * @param a
     *  creation spin orbital index
     * @param i
     *  annihilation spin orbital index
     * @return
     *  one-electron coefficient of H given by the 2D array T
     */
    virtual defs::ham_t get_coeff_1100(size_t a, size_t i) const {return 0;}
    /**
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

    defs::ham_t get_element(const field::FrmOnv &ket, const conn::FrmOnv &conn) const;

    size_t nci() const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    virtual void log_data() const;

    virtual bool enabled() const {
        return true;
    }

    bool disabled() const {
        return !enabled();
    }
};

struct NullFrmHam : FrmHam {
    NullFrmHam() : FrmHam(0, 0, true, {}){}

    bool enabled() const override {
        return false;
    }
};

/**
 * fermion sites may not be doubly occupied or unoccupied in spin systems
 */
struct SpinModelFrmHam : FrmHam {
    SpinModelFrmHam(size_t nelec, size_t nsite, int ms2_restrict): FrmHam(nelec, nsite, ms2_restrict){}
};

#endif //M7_FRMHAM_H
