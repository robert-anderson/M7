//
// Created by anderson on 12/9/21.
//

#ifndef M7_GENERALLADDERHAM_H
#define M7_GENERALLADDERHAM_H

#include <src/core/io/EbdumpFileReader.h>
#include <src/core/integrals/FrmBosCoupledCoeffs.h>
#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

#if 0
struct GeneralLadderHam {

    const size_t m_nboson_max;
    const BasisDims m_bd;
    /**
     * coefficients for "coupled" ranksigs 1110, 1101. contributing exsigs are either:
     *  "density coupled" (0010, 0001), or
     *  "hopping coupled" (1110, 1101)
     */
    FrmBosCoupledCoeffs m_v;
    /**
     * coefficients for "uncoupled" ranksigs 0010, 0001. only contributing exsigs are
     *  "uncoupled" (0010, 0001)
     *
     * density-coupled and uncoupled excitations have the same exsig, collectively they will be called "pure" bosonic
     * excitations / de-excitations, as opposed to the fermion-coupled "hopping" exsigs
     */
    std::vector<defs::ham_t> m_v_unc;

    ham_data::TermContribs m_contribs_0010;
    ham_data::TermContribs m_contribs_0001;
    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    GeneralLadderHam(const EbdumpHeader& header, size_t nboson_max);

    GeneralLadderHam(const std::string& fname, size_t nboson_max) :
            LadderHam(EbdumpHeader(fname), nboson_max){}

    defs::ham_t get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;
};

#endif //M7_GENERALLADDERHAM_H
#endif //M7_GENERALLADDERHAM_H