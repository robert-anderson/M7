//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminants.h"

void OccupiedUpdater::operator()(const field::FrmOnv &onv, defs::inds &inds) {
    DEBUG_ASSERT_LE(onv.nbit(), inds.capacity(), "occupied updater inds not large enough for ONV");
    inds.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void VacantUpdater::operator()(const field::FrmOnv &onv, defs::inds &inds) {
    DEBUG_ASSERT_LE(onv.nbit(), inds.capacity(), "vacant updater inds not large enough for ONV");
    inds.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}


void NdOccupiedUpdater::operator()(const field::FrmOnv &onv, const defs::inds& map,
        defs::inds& flat_inds, std::vector<defs::inds> &nd_inds) {
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}


void NdVacantUpdater::operator()(const field::FrmOnv &onv, const defs::inds& map,
        defs::inds& flat_inds, std::vector<defs::inds> &nd_inds){
    flat_inds.clear();
    for (auto& v :nd_inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            nd_inds[map[ibit]].push_back(ibit);
            flat_inds.push_back(ibit);
        }
    }
}