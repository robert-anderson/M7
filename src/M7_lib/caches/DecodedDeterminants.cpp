//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminants.h"

#if 0
void OccupiedUpdater::operator()(const field::FrmOnv &onv, uintv_t &uintv_t) {
    DEBUG_ASSERT_LE(onv.nbit(), uintv_t.capacity(), "occupied updater uintv_t not large enough for ONV");
    uintv_t.clear();
    for (uint_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) uintv_t.push_back(bit::next_setbit(work) + idataword * nbit_word);
    }
}

void VacantUpdater::operator()(const field::FrmOnv &onv, uintv_t &uintv_t) {
    DEBUG_ASSERT_LE(onv.nbit(), uintv_t.capacity(), "vacant updater uintv_t not large enough for ONV");
    uintv_t.clear();
    for (uint_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) uintv_t.push_back(bit::next_setbit(work) + idataword * nbit_word);
    }
}


void NdOccupiedUpdater::operator()(const field::FrmOnv &onv, const uintv_t& map,
        uintv_t& flat_inds, v_t<uintv_t> &nd_inds) {
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (uint_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) {
            auto ibit = bit::next_setbit(work) + idataword * nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}


void NdVacantUpdater::operator()(const field::FrmOnv &onv, const uintv_t& map,
        uintv_t& flat_inds, v_t<uintv_t> &nd_inds){
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (uint_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) {
            auto ibit = bit::next_setbit(work) + idataword * nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}
#endif