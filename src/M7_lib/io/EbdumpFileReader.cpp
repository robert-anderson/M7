//
// Created by rja on 22/08/2021.
//

#include "EbdumpFileReader.h"


EbdumpHeader::EbdumpHeader(const std::string &fname) :
        FortranNamelistReader(fname), m_nmode(read_int("NMODE")), m_nsite(read_int("NORB")), m_uhf(read_bool("UHF")) {}


EbdumpFileReader::EbdumpFileReader(const std::string &fname) :
        HamiltonianFileReader(fname, 3),
        m_header(fname), m_norb_distinct((m_header.m_uhf ? 2 : 1)*m_header.m_nsite){}

size_t EbdumpFileReader::ranksig(const defs::inds &inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
    switch (nset_ind(inds)) {
        case 1ul:
            DEBUG_ASSERT_NE(inds[0], ~0ul, "set index should be the first one: the boson mode index");
            return exsig_utils::ex_0010;
        case 3ul:
            return exsig_utils::ex_1110;
        default:
            return ~0ul;
    }
}

size_t EbdumpFileReader::exsig(const defs::inds &inds, const size_t &ranksig) const {
    DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
    if (inds[1] == inds[2]) return exsig_utils::ex_0010;
    else return exsig_utils::ex_1110;
}

bool EbdumpFileReader::inds_in_range(const defs::inds &inds) const {
    return inds[0] <= m_header.m_nmode && inds[1] <= m_norb_distinct && inds[2] <= m_norb_distinct;
}
