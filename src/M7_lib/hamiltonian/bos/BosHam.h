//
// Created by Robert J. Anderson on 26/07/2021.
//

#ifndef M7_BOSHAM_H
#define M7_BOSHAM_H

#include "M7_lib/integrals/BosonCoeffs_1.h"
#include "M7_lib/integrals/BosonCoeffs_2.h"
#include "M7_lib/io/BosdumpFileReader.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/parallel/SharedArray.h"
#include "M7_lib/field/Fields.h"

#include "M7_lib/hamiltonian/HamiltonianData.h"
#include "M7_lib/hamiltonian/HamOpTerm.h"

struct BosHam : HamOpTerm {
    /**
     * a convenient pair of references to the relevant Hamiltonian section and the Basis configuration section
     */
    typedef HamOpTerm::OptPair<conf::BosHam> opt_pair_t;
    /**
     * properties of the many-body basis
     */
    const sys::bos::Basis m_basis;
    /**
     * term contributions. the four digit represent the rank signature.
     * e.g. hamiltonian term of rank 0001 can only take contributions from excitations of exsig 0001.
     * on the other hand, terms of rank 0022 can take contributions from exsigs 0000 (diagonals), 0011, and 0022
     * these objects keep track of which of these exsigs are non-zero and which may contribute to matrix elements
     */
    ham::TermContribs m_contribs_0010;
    ham::TermContribs m_contribs_0001;
    ham::TermContribs m_contribs_0011;
    ham::TermContribs m_contribs_0022;

    BosHam(const sys::bos::Basis& basis):
            m_basis(basis),
            m_contribs_0010(exsig::ex_0010), m_contribs_0001(exsig::ex_0001),
            m_contribs_0011(exsig::ex_0011), m_contribs_0022(exsig::ex_0022) {}

	virtual ~BosHam(){}

    virtual ham_t get_coeff_0011(uint_t /*i*/, uint_t /*j*/) const {return 0;}
    virtual ham_t get_coeff_0022(uint_t /*i*/, uint_t /*j*/,
                                       uint_t /*k*/, uint_t /*l*/) const {return 0;}

    virtual ham_t get_element_0000(const field::BosOnv& /*onv*/) const {return 0;}
    virtual ham_t get_element_0011(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const {return 0;}
    virtual ham_t get_element_0022(const field::BosOnv& /*onv*/, const conn::BosOnv& /*conn*/) const {return 0;}

    ham_t get_element(const field::BosOnv& onv) const {
        return get_element_0000(onv);
    }

    ham_comp_t get_energy(const field::BosOnv& onv) const {
        return arith::real(get_element(onv));
    }

    ham_t get_element(const field::BosOnv& src, const conn::BosOnv& conn) const {
        switch (conn.size()) {
            case 0: return get_element_0000(src);
            case 2: return get_element_0011(src, conn);
            case 4: return get_element_0022(src, conn);
        }
        return 0;
    }

    virtual void log_data() const;

    virtual uint_t default_nboson() const {
        return 0ul;
    }

};

/**
 * boson hamiltonian which can be defined in a non-zero number of modes, but with no non-zero term coefficients
 */
struct NullBosHam : BosHam, NullOpTerm {
    NullBosHam() : BosHam(0ul){}
};

#endif //M7_BOSHAM_H