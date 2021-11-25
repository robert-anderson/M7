//
// Created by rja on 26/07/2021.
//

#ifndef M7_BOSONHAMILTONIAN_H
#define M7_BOSONHAMILTONIAN_H

#include <src/core/integrals/BosonCoeffs_1.h>
#include <src/core/integrals/BosonCoeffs_2.h>
#include <src/core/io/BosdumpFileReader.h>
#include "src/core/connection/Connections.h"
#include "src/core/parallel/SharedArray.h"
#include "src/core/field/Fields.h"
#include "HamiltonianData.h"

struct BosonHamiltonian {
    const size_t m_nmode, m_nboson, m_nboson_max;
    BosonCoeffs_1 m_coeffs_1;
    BosonCoeffs_2 m_coeffs_2;
    ham_data::TermContribs m_contribs_0011;
    ham_data::TermContribs m_contribs_0022;

    BosonHamiltonian(const std::string& fname, size_t nboson_max);

    defs::ham_t get_element(const field::BosOnv &onv) const;

    defs::ham_comp_t get_energy(const field::BosOnv &onv) const;

    defs::ham_t get_element(const field::BosOnv &onv, const conn::BosOnv& conn) const;

    size_t nci() const;

    void log_data() const;

private:
    static size_t read_nmode(const std::string& fname){
        if (!FileReader::exists(fname)) return 0ul;
        return BosdumpFileReader(fname).m_nmode;
    }
    static size_t read_nboson(const std::string& fname){
        if (!FileReader::exists(fname)) return 0ul;
        return BosdumpFileReader(fname).m_nboson;
    }
};


#endif //M7_BOSONHAMILTONIAN_H
