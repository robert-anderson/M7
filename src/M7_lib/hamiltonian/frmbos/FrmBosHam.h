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

struct FrmBosHam {

    const BasisData m_bd;
    const size_t m_nboson_max;

    ham_data::TermContribs m_contribs_0010;
    ham_data::TermContribs m_contribs_0001;
    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    FrmBosHam(const BasisData &bd, size_t nboson_max);
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
    NullLadderHam() : FrmBosHam({0, 0}, 0){}

    bool enabled() const override {
        return false;
    }
};


#endif //M7_FRMBOSHAM_H
