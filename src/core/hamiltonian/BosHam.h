//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSHAM_H
#define M7_BOSHAM_H

#include <src/core/integrals/BosonCoeffs_1.h>
#include <src/core/integrals/BosonCoeffs_2.h>
#include <src/core/io/BosdumpFileReader.h>
#include "src/core/connection/Connections.h"
#include "src/core/parallel/SharedArray.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct BosHam {
    const size_t m_nmode, m_nboson;
    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    BosHam(size_t nmode, size_t nboson);

    virtual defs::ham_t get_coeff_0011(const size_t& i, const size_t& j) const {return 0;}
    virtual defs::ham_t get_coeff_0022(const size_t& i, const size_t& j,
                                       const size_t& k, const size_t& l) const {return 0;}

    virtual defs::ham_t get_element_0000(const field::BosOnv &onv) const {return 0;}
    virtual defs::ham_t get_element_0011(const field::BosOnv &onv, const conn::BosOnv& conn) const {return 0;}
    virtual defs::ham_t get_element_0022(const field::BosOnv &onv, const conn::BosOnv& conn) const {return 0;}

    defs::ham_t get_element(const field::BosOnv &onv) const {
        return get_element_0000(onv);
    }

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const {
        return consts::real(get_element(onv));
    }

    defs::ham_t get_element(const field::BosOnv &src, const conn::BosOnv& conn) const {
        switch (conn.size()) {
            case 0: return get_element_0000(src);
            case 2: return get_element_0011(src, conn);
            case 4: return get_element_0022(src, conn);
        }
        return 0;
    }

    size_t nci() const;

    void log_data() const;

};


#endif //M7_BOSHAM_H