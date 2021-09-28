//
// Created by rja on 05/11/2020.
//

#ifndef M7_LADDERHAMILTONIAN_H
#define M7_LADDERHAMILTONIAN_H

#include <src/core/io/EbdumpFileReader.h>
#include <src/core/integrals/FrmBosCoupledCoeffs.h>
#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct LadderHamiltonian {

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

    LadderHamiltonian(const std::string& fname, size_t nboson_max);

    defs::ham_t get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const;

    bool is_holstein() const;

    bool constant_uncoupled() const;

    bool is_zpm_half_filled() const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;

private:
    static BasisDims read_bd(const std::string& fname);

    static bool read_spin_resolved(const std::string& fname);
};

#endif //M7_LADDERHAMILTONIAN_H