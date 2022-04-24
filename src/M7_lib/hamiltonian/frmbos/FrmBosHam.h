//
// Created by rja on 05/11/2020.
//

#ifndef M7_FRMBOSHAM_H
#define M7_FRMBOSHAM_H

#include "M7_lib/io/EbdumpFileReader.h"
#include "M7_lib/integrals/FrmBosCoupledCoeffs.h"
#include "M7_lib/connection/Connections.h"
#include "M7_lib/field/Fields.h"

#include "M7_lib/hamiltonian/HamiltonianData.h"
#include "M7_lib/hamiltonian/HamOpTerm.h"
#include "M7_lib/hamiltonian/frm/FrmHam.h"
#include "M7_lib/hamiltonian/bos/BosHam.h"

/**
 * base class for all Hamiltonians expressed in terms of products of fermionic and bosonic second-quantised operators
 */
struct FrmBosHam : HamOpTerm {

    const HilbertSpace m_hs;

    ham_data::TermContribs m_contribs_0010;
    ham_data::TermContribs m_contribs_0001;
    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    /**
     * @param bd
     *  basis data determined by the information available to the derived class ctors. its compatibility with any
     *  defined fermion or boson ham is checked here by referencing the frm and bos arguments in contrast, the
     *  fermion-boson product term does not determine anything about the Hilbert space, and so these attributes are
     *  brought in directly from the fermion and boson parts.
     * @param frm
     *  fermionic part of H, from which the fermionic Hilbert space attributes are copied
     * @param bos
     *  bosonic part of H, from which the bosonic Hilbert space attributes are copied
     */
    FrmBosHam(const HilbertSpace& hs, const FrmHam& frm, const BosHam& bos);
    /*
     * all Hilbert space info is taken directly from the FrmHam and BosHam
     */
    FrmBosHam(const FrmHam& frm, const BosHam& bos);

    virtual ~FrmBosHam(){}

    virtual defs::ham_t get_coeff_0010(size_t imode) const {return 0;}
    virtual defs::ham_t get_coeff_0001(size_t imode) const {return 0;}
    virtual defs::ham_t get_coeff_1110(size_t imode, size_t i, size_t j) const {return 0;}
    virtual defs::ham_t get_coeff_1101(size_t imode, size_t i, size_t j) const {return 0;}

    virtual defs::ham_t get_element_0010(const field::BosOnv& onv, const conn::BosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_0001(const field::BosOnv& onv, const conn::BosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_0010(const field::FrmBosOnv& onv, const conn::FrmBosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_0001(const field::FrmBosOnv& onv, const conn::FrmBosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_1110(const field::FrmBosOnv& onv, const conn::FrmBosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_1101(const field::FrmBosOnv& onv, const conn::FrmBosOnv& conn) const {return 0;}


    defs::ham_t get_element(const field::BosOnv &onv, const conn::BosOnv &conn) const {
        switch (conn.exsig()) {
            case exsig_utils::ex_0001: return get_element_0001(onv, conn);
            case exsig_utils::ex_0010: return get_element_0010(onv, conn);
        }
        return 0.0;
    }

    defs::ham_t get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        switch (conn.exsig()) {
            case exsig_utils::ex_0001: return get_element_0001(onv, conn);
            case exsig_utils::ex_0010: return get_element_0010(onv, conn);
            case exsig_utils::ex_1101: return get_element_1101(onv, conn);
            case exsig_utils::ex_1110: return get_element_1110(onv, conn);
        }
        return 0.0;
    }

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;

    virtual bool enabled() const {
        return true;
    }

    bool disabled() const {
        return !enabled();
    }
};

struct NullLadderHam: FrmBosHam {
    NullLadderHam() : FrmBosHam({}, NullFrmHam(), NullBosHam()){}

    bool enabled() const override {
        return false;
    }
};


#endif //M7_FRMBOSHAM_H
