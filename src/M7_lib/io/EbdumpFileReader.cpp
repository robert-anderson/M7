//
// Created by Robert J. Anderson on 22/08/2021.
//

#include "EbdumpFileReader.h"
#include "M7_lib/util/Exsig.h"


EbdumpFileReader::EbdumpFileReader(const EbdumpInfo& info) :
        HamTextFileReader(info.m_fname, 3), m_info(info),
        m_norb_distinct((m_info.m_uhf ? 2 : 1)*m_info.m_nsite) {}

uint_t EbdumpFileReader::ranksig(const uintv_t& inds) const {
    DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
    switch (nset_ind(inds)) {
        case 1ul:
            DEBUG_ASSERT_NE(inds[0], ~0ul, "set index should be the first one: the boson mode index");
            return exsig::ex_0010;
        case 3ul:
            return exsig::ex_1110;
        default:
            return ~0ul;
    }
}

uint_t EbdumpFileReader::exsig(const uintv_t& inds, uint_t) const {
    DEBUG_ASSERT_EQ(inds.size(), 3ul, "incorrect maximum number of SQ operator indices");
    if (inds[1] == inds[2]) return exsig::ex_0010;
    else return exsig::ex_1110;
}

bool EbdumpFileReader::inds_in_range(const uintv_t& inds) const {
    return inds[0] <= m_info.m_nmode && inds[1] <= m_norb_distinct && inds[2] <= m_norb_distinct;
}

void EbdumpFileReader::convert_inds(uintv_t& inds) {
    if (!m_info.m_spin_resolved || m_info.m_spin_major) return;
    auto fn = [&inds, this](uint_t i){
        auto& ind = inds[i];
        ind = (ind == ~0ul) ? ~0ul : (ind / 2 + ((ind & 1ul) ? m_info.m_nsite : 0));
    };
    // this conversion only applies to the fermion indices
    fn(1);
    fn(2);
}

bool EbdumpFileReader::next(uintv_t& inds, ham_t& v) {
    if (!HamTextFileReader::next(inds, v)) return false;
    convert_inds(inds);
    return true;
}
