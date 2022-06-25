//
// Created by Robert J. Anderson on 19/08/2021.
//

#ifndef M7_BOSDUMPFILEREADER_H
#define M7_BOSDUMPFILEREADER_H

#include "HamTextFileReader.h"
#include "FortranNamelistReader.h"

struct BosdumpHeader : FortranNamelistReader {
    const uint_t m_nmode, m_nboson;
    BosdumpHeader(const std::string& fname);
};

struct BosdumpFileReader : HamTextFileReader {
    BosdumpHeader m_header;

    BosdumpFileReader(const std::string &fname);

    uint_t ranksig(const uintv_t &inds) const override;

    uint_t exsig(const uintv_t &inds, uint_t ranksig) const override;

    bool inds_in_range(const uintv_t &inds) const override;
};


#endif //M7_BOSDUMPFILEREADER_H
