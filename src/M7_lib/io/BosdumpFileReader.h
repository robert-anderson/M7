//
// Created by Robert J. Anderson on 19/08/2021.
//

#ifndef M7_BOSDUMPFILEREADER_H
#define M7_BOSDUMPFILEREADER_H

#include "HamiltonianFileReader.h"
#include "FortranNamelistReader.h"

struct BosdumpHeader : FortranNamelistReader {
    const size_t m_nmode, m_nboson;
    BosdumpHeader(const std::string& fname);
};

struct BosdumpFileReader : HamiltonianFileReader {
    BosdumpHeader m_header;

    BosdumpFileReader(const std::string &fname);

    size_t ranksig(const defs::ivec_t &inds) const override;

    size_t exsig(const defs::ivec_t &inds, size_t ranksig) const override;

    bool inds_in_range(const defs::ivec_t &inds) const override;
};


#endif //M7_BOSDUMPFILEREADER_H
