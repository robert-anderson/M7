//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "BosdumpFileReader.h"

BosdumpHeader::BosdumpHeader(const str_t &fname) :
        FortranNamelistReader(fname), m_nmode(read_int("NMODE")), m_nboson(read_int("NBOSON")) {}

BosdumpFileReader::BosdumpFileReader(const str_t &fname) :
        HamTextFileReader(fname, 4), m_header(BosdumpHeader(fname)) {}

OpSig BosdumpFileReader::ranksig(const uintv_t &inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    return inds[2] == ~0ul ? opsig::c_0011 : opsig::c_0022;
}

OpSig BosdumpFileReader::exsig(const uintv_t &inds, OpSig ranksig) const {
    DEBUG_ASSERT_EQ(inds.size(), 4ul, "incorrect maximum number of SQ operator indices");
    switch (ranksig) {
        case opsig::c_zero:
            return opsig::c_zero;
        case opsig::c_0011:
            return inds[0] == inds[1] ? opsig::c_zero : opsig::c_0011;
        case opsig::c_0022:
            return inds[0] == inds[1] ?
                   (inds[2] == inds[3] ? opsig::c_zero : opsig::c_0011) :
                   (inds[2] == inds[3] ? opsig::c_0011 : opsig::c_0022);
        default:
            return opsig::c_invalid;
    }
}

bool BosdumpFileReader::inds_in_range(const uintv_t &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](uint_t i) { return i <= m_header.m_nmode; });
}