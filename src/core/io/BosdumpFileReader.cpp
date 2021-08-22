//
// Created by rja on 19/08/2021.
//

#include "BosdumpFileReader.h"

BosdumpFileReader::BosdumpFileReader(const std::string &fname) : HamiltonianFileReader(fname, 3, false){
    REQUIRE_FALSE_ALL(m_spin_resolved, "spin resolved boson dumps are invalid");
}

size_t BosdumpFileReader::ranksig(const defs::inds &inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 2ul, "incorrect maximum number of SQ operator indices");
    return conn_utils::encode_exsig(0, 0, 1, 1);
}

size_t BosdumpFileReader::exsig(const defs::inds &inds, const size_t ranksig) const {
    DEBUG_ASSERT_EQ(inds.size(), 2ul, "incorrect maximum number of SQ operator indices");
    size_t n = inds[0]!=inds[1];
    return conn_utils::encode_exsig(0, 0, n, n);
}
