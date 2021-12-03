//
// Created by rja on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "HamiltonianFileReader.h"
#include "FortranNamelistReader.h"

struct EbdumpHeader : FortranNamelistReader {
    const size_t m_nmode, m_nsite;
    const bool m_uhf;
    EbdumpHeader(const std::string& fname);
};

struct EbdumpFileReader : HamiltonianFileReader {
    const EbdumpHeader m_header;
    const size_t m_norb_distinct;

    EbdumpFileReader(const std::string &fname);

    size_t ranksig(const defs::inds &inds) const override;

    size_t exsig(const defs::inds &inds, const size_t &ranksig) const override;

    bool inds_in_range(const defs::inds &inds) const override;
};


#endif //M7_EBDUMPFILEREADER_H
