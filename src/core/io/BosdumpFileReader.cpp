//
// Created by rja on 19/08/2021.
//

#include "BosdumpFileReader.h"


BosdumpFileReader::BosdumpFileReader(const std::string &fname) : HamiltonianFileReader(fname, 4), m_header(BosdumpHeader(fname)){}


size_t BosdumpFileReader::ranksig(const defs::inds &inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    return inds[2]==~0ul ? exsig_utils::ex_0011 : exsig_utils::ex_0022;
}

size_t BosdumpFileReader::exsig(const defs::inds &inds, const size_t &ranksig) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    switch (ranksig) {
        case 0ul:
            return 0ul;
        case exsig_utils::ex_0011:
            return inds[0] == inds[1] ? 0ul : exsig_utils::ex_0011;
        case exsig_utils::ex_0022:
            return inds[0] == inds[1] ?
                   (inds[2] == inds[3] ? 0ul : exsig_utils::ex_0011) :
                   (inds[2] == inds[3] ? exsig_utils::ex_0011 : exsig_utils::ex_0022);
        default:
            return ~0ul;
    }
}

bool BosdumpFileReader::inds_in_range(const defs::inds &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](size_t i){return i<=m_header.m_nmode;});
}