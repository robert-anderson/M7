//
// Created by rja on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "HamiltonianFileReader.h"

struct EbdumpFileReader : HamiltonianFileReader {
    const size_t m_nmode;

    EbdumpFileReader(const std::string &fname): HamiltonianFileReader(fname, 3, false),
        m_nmode(read_header_int(fname, "NMODE")){
        REQUIRE_FALSE_ALL(m_spin_resolved, "spin resolved electron-boson dumps are not currently supported");
    }

    size_t ranksig(const defs::inds &inds) const override {
        DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
        switch (nset_ind(inds)) {
            case 1ul:
                DEBUG_ASSERT_NE(inds[0], ~0ul, "set index should be the first one: the boson mode index");
                return exsig_utils::ex_0010;
            case 3ul: return exsig_utils::ex_1110;
            default: return ~0ul;
        }
    }

    size_t exsig(const defs::inds &inds, const size_t& ranksig) const override {
        DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
        if (inds[1]==inds[2]) return exsig_utils::ex_0010;
        else return exsig_utils::ex_1110;
    }

};


#endif //M7_EBDUMPFILEREADER_H
