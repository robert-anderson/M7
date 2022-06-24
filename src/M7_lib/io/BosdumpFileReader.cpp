//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "BosdumpFileReader.h"
#include <M7_lib/util/Exsig.h>

BosdumpHeader::BosdumpHeader(const std::string &fname) :
        FortranNamelistReader(fname), m_nmode(read_int("NMODE")), m_nboson(read_int("NBOSON")) {}

BosdumpFileReader::BosdumpFileReader(const std::string &fname) :
        HamiltonianFileReader(fname, 4), m_header(BosdumpHeader(fname)) {}

uint_t BosdumpFileReader::ranksig(const defs::uintv_t &inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    return inds[2] == ~0ul ? exsig::ex_0011 : exsig::ex_0022;
}

uint_t BosdumpFileReader::exsig(const defs::uintv_t &inds, uint_t ranksig) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    switch (ranksig) {
        case 0ul:
            return 0ul;
        case exsig::ex_0011:
            return inds[0] == inds[1] ? 0ul : exsig::ex_0011;
        case exsig::ex_0022:
            return inds[0] == inds[1] ?
                   (inds[2] == inds[3] ? 0ul : exsig::ex_0011) :
                   (inds[2] == inds[3] ? exsig::ex_0011 : exsig::ex_0022);
        default:
            return ~0ul;
    }
}

bool BosdumpFileReader::inds_in_range(const defs::uintv_t &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](uint_t i) { return i <= m_header.m_nmode; });
}