//
// Created by rja on 27/02/2020.
//

#ifndef M7_FRMHAM_H
#define M7_FRMHAM_H

#include <cstddef>

#include <M7_lib/caches/DecodedDeterminants.h>
#include <M7_lib/connection/Connections.h>
#include <M7_lib/io/Options.h>
#include <M7_lib/config/FciqmcConfig.h>
#include <M7_lib/integrals/Integrals_1e.h>
#include <M7_lib/integrals/Integrals_2e.h>
#include <M7_lib/table/BufferedFields.h>

#include "HamiltonianData.h"

/**
 * All interactions between the fermionic parts of MBFs are described in this class.
 */
struct FrmHam {

    const size_t m_nelec;
    const size_t m_nsite;
    const int m_ms2_restrict;
    const bool m_complex_valued = false;

    AbelianGroupMap m_point_group_map{PointGroup(), defs::inds(m_nsite, 0ul)};

    defs::ham_t m_e_core = 0.0;

    ham_data::TermContribs m_contribs_1100;
    ham_data::TermContribs m_contribs_2200;
    ham_data::KramersAttributes m_kramers_attrs;

    FrmHam() = delete;

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict,
           bool complex_valued, defs::inds site_irreps);

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict, bool complex_valued);

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict, defs::inds site_irreps);

    FrmHam(size_t nelec, size_t nsite, int ms2_restrict);

    FrmHam(const FrmHam& other)
        : FrmHam(other.m_nelec, other.m_nsite, other.m_ms2_restrict,
                 other.m_complex_valued, other.m_point_group_map.m_site_irreps){}

    FrmHam& operator=(const FrmHam& other){return *this;}

	virtual ~FrmHam(){}

    virtual defs::ham_t get_coeff_1100(size_t i, size_t j) const {return 0;}
    virtual defs::ham_t get_coeff_2200(size_t i, size_t j,
                                       size_t k, size_t l) const {return 0;}

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

    virtual bool enabled() const {
        return true;
    }

    bool disabled() const {
        return !enabled();
    }
};

struct NullFrmHam : FrmHam {
    NullFrmHam() : FrmHam(0, 0, true, false, {}){}

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
