//
// Created by rja on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "HamiltonianFileReader.h"

struct EbdumpFileReader : HamiltonianFileReader {
    const size_t m_nmode;

    EbdumpFileReader(const std::string &fname);

    size_t ranksig(const defs::inds &inds) const override;

    size_t exsig(const defs::inds &inds, const size_t &ranksig) const override;

    bool inds_in_range(const defs::inds &inds) const override;
};


#endif //M7_EBDUMPFILEREADER_H
