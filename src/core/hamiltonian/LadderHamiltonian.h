//
// Created by rja on 05/11/2020.
//

#ifndef M7_LADDERHAMILTONIAN_H
#define M7_LADDERHAMILTONIAN_H

#include <src/core/io/EbdumpFileReader.h>
#include "src/core/connection/Connections.h"
#include "src/core/field/Fields.h"
#include "src/core/parallel/SharedArray.h"
#include "HamiltonianData.h"

class LadderHamiltonian {
    size_t index(const size_t &n, const size_t &p, const size_t &q) const {
        return n * m_nmode2 + p * m_nmode + q;
    }

public:

    const size_t m_nboson_max, m_nmode;
    const size_t m_nmode2;
    SharedArray<defs::ham_t> m_v;

    ham_data::TermContribs m_contribs_1110;
    ham_data::TermContribs m_contribs_1101;

    LadderHamiltonian(size_t nmode, size_t nboson_max, std::string fname);

    defs::ham_t v(const size_t &n, const size_t &p, const size_t &q) const;

    defs::ham_t get_element(const field::FrmBosOnv &onv, const conn::FrmBosOnv &conn) const;

    bool is_holstein() const {
        return !m_contribs_1101.is_nonzero(exsig_utils::ex_1101);
    }

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;

};

#endif //M7_LADDERHAMILTONIAN_H