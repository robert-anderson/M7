//
// Created by Robert J. Anderson on 05/11/2020.
//

#ifndef M7_FRMBOSHAM_H
#define M7_FRMBOSHAM_H

#include "M7_lib/io/EbdumpFileReader.h"
#include "M7_lib/integrals/FrmBosCoupledCoeffs.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/field/Fields.h"

#include "M7_lib/hamiltonian/HamiltonianData.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"
#include "M7_lib/hamiltonian/bos/BosHam.h"

/**
 * base class for all Hamiltonians expressed in terms of products of fermionic and bosonic second-quantised operators
 */
struct FrmBosHam : HamOpTerm {
    /**
     * a convenient bundle of references to the relevant initialization information
     */
    typedef HamOpTerm::InitOpts<conf::FrmBosHam> init_opts_t;

    const sys::Basis m_basis;

    /**
     * term contributions. the four digit represent the rank signature.
     * e.g. hamiltonian term of rank 1101 can take contributions from excitations of exsig 0001 and 1101.
     * these objects keep track of which of these exsigs are non-zero and which may contribute to matrix elements
     */
    ham::TermContribs m_contribs_1110;
    ham::TermContribs m_contribs_1101;

    FrmBosHam(sys::Basis basis);

    virtual ~FrmBosHam(){}

    virtual ham_t get_coeff_1110(uint_t /*imode*/, uint_t /*i*/, uint_t /*j*/) const {return 0;}
    virtual ham_t get_coeff_1101(uint_t /*imode*/, uint_t /*i*/, uint_t /*j*/) const {return 0;}

    virtual ham_t get_element_0010(const field::FrmBosOnv& /*onv*/,
                                         const conn::FrmBosOnv& /*conn*/) const {return 0;}
    virtual ham_t get_element_0001(const field::FrmBosOnv& /*onv*/,
                                         const conn::FrmBosOnv& /*conn*/) const {return 0;}
    virtual ham_t get_element_1110(const field::FrmBosOnv& /*onv*/,
                                         const conn::FrmBosOnv& /*conn*/) const {return 0;}
    virtual ham_t get_element_1101(const field::FrmBosOnv& /*onv*/,
                                         const conn::FrmBosOnv& /*conn*/) const {return 0;}


    ham_t get_element(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const {
        return 0.0;
    }

    ham_t get_element(const field::FrmBosOnv& onv, const conn::FrmBosOnv& conn) const {
        switch (conn.exsig()) {
            case exsig::ex_0001: return get_element_0001(onv, conn);
            case exsig::ex_0010: return get_element_0010(onv, conn);
            case exsig::ex_1101: return get_element_1101(onv, conn);
            case exsig::ex_1110: return get_element_1110(onv, conn);
        }
        return 0.0;
    }

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;

    virtual uint_t default_nelec() const {
        // assume 1/2-filling
        return m_basis.m_frm.m_nsite;
    }

    virtual int default_ms2_value() const {
        return sys::frm::c_undefined_ms2;
    }

    virtual uint_t default_nboson() const {
        return 0ul;
    }
};

/**
 * fermion-boson hamiltonian which can be defined in non-zero numbers of sites and modes, but with no non-zero term
 * coefficients
 */
struct NullFrmBosHam : FrmBosHam, NullOpTerm {
    NullFrmBosHam() : FrmBosHam({0ul, 0ul}){
        DEBUG_ASSERT_FALSE(m_contribs_1101.is_nonzero(exsig::ex_0001) || m_contribs_1101.is_nonzero(exsig::ex_1101),
                           "No exsigs should contribute to null FrmBosHam");
        DEBUG_ASSERT_FALSE(m_contribs_1110.is_nonzero(exsig::ex_0010) || m_contribs_1110.is_nonzero(exsig::ex_1110),
                           "No exsigs should contribute to null FrmBosHam");
    }
};

#endif //M7_FRMBOSHAM_H