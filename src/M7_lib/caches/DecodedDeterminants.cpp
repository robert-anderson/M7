//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminants.h"

#if 0
void OccupiedUpdater::operator()(const field::FrmOnv &onv, defs::ivec_t &ivec_t) {
    DEBUG_ASSERT_LE(onv.nbit(), ivec_t.capacity(), "occupied updater ivec_t not large enough for ONV");
    ivec_t.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) ivec_t.push_back(bit::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void VacantUpdater::operator()(const field::FrmOnv &onv, defs::ivec_t &ivec_t) {
    DEBUG_ASSERT_LE(onv.nbit(), ivec_t.capacity(), "vacant updater ivec_t not large enough for ONV");
    ivec_t.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) ivec_t.push_back(bit::next_setbit(work) + idataword * defs::nbit_word);
    }
}


void NdOccupiedUpdater::operator()(const field::FrmOnv &onv, const defs::ivec_t& map,
        defs::ivec_t& flat_inds, std::vector<defs::ivec_t> &nd_inds) {
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) {
            auto ibit = bit::next_setbit(work) + idataword * defs::nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}


void NdVacantUpdater::operator()(const field::FrmOnv &onv, const defs::ivec_t& map,
        defs::ivec_t& flat_inds, std::vector<defs::ivec_t> &nd_inds){
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) {
            auto ibit = bit::next_setbit(work) + idataword * defs::nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}
#endif