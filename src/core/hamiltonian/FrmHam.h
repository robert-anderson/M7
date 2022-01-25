//
// Created by rja on 27/02/2020.
//

#ifndef M7_FRMHAM_H
#define M7_FRMHAM_H

#include <cstddef>
#include <src/core/basis/DecodedDeterminants.h>
#include <src/core/connection/Connections.h>
#include <src/core/io/Options.h>
#include <src/core/config/FciqmcConfig.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "src/core/table/BufferedFields.h"
#include "HamiltonianData.h"

/**
 * All interactions between the fermionic parts of MBFs are described in this class.
 */
struct FrmHam {

    const size_t m_nelec;
    const size_t m_nsite;
    const int m_ms2_restrict;
    const bool m_complex_valued;

    AbelianGroupMap m_point_group_map;

    defs::ham_t m_e_core = 0.0;

    ham_data::TermContribs m_contribs_1100;
    ham_data::TermContribs m_contribs_2200;
    ham_data::KramersAttributes m_kramers_attrs;

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict,
                       bool complex_valued=false, defs::inds site_irreps = {});
	virtual ~FrmHam(){}

    virtual defs::ham_t get_coeff_1100(const size_t& i, const size_t& j) const {return 0;}
    virtual defs::ham_t get_coeff_2200(const size_t& i, const size_t& j,
                                       const size_t& k, const size_t& l) const {return 0;}

    virtual defs::ham_t get_element_0000(const field::FrmOnv &onv) const {return 0;}
    virtual defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {return 0;}
    virtual defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {return 0;}
    //virtual defs::ham_t get_element_3300(const field::FrmOnv &onv, const conn::FrmOnv &conn) const = 0;

    defs::ham_t get_element(const field::FrmOnv &onv) const;

    defs::ham_comp_t get_energy(const field::FrmOnv &onv) const;

    defs::ham_t get_element(const field::FrmOnv &onv, const conn::FrmOnv &conn) const;

    size_t nci() const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    virtual void log_data() const;
};

#endif //M7_FRMHAM_H
