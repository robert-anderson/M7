//
// Created by rja on 05/11/2020.
//

#ifndef M7_LADDERHAM_H
#define M7_LADDERHAM_H

#include <src/core/io/EbdumpFileReader.h>
#include <src/core/integrals/FrmBosCoupledCoeffs.h>
#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct LadderHam {

    const BasisDims m_bd;
    const size_t m_nboson_max;

    ham_data::TermContribs m_contribs_0010;
    ham_data::TermContribs m_contribs_0001;
    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    LadderHam(const BasisDims &bd, size_t nboson_max);
	virtual ~LadderHam(){}

    virtual defs::ham_t get_coeff_0010(const size_t& imode) const {return 0;}
    virtual defs::ham_t get_coeff_0001(const size_t& imode) const {return 0;}
    virtual defs::ham_t get_coeff_1110(const size_t &imode, const size_t &j, const size_t &i) const {return 0;}
    virtual defs::ham_t get_coeff_1101(const size_t &imode, const size_t &j, const size_t &i) const {return 0;}

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
};

#endif //M7_LADDERHAM_H
