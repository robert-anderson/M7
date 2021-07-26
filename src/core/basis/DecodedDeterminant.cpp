//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

void OccupiedUpdater::operator()(const fields::FrmOnv &onv, defs::inds &inds) {
    DEBUG_ASSERT_LE(onv.nbit(), inds.capacity(), "occupied updater inds not large enough for ONV");
    inds.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void VacantUpdater::operator()(const fields::FrmOnv &onv, defs::inds &inds) {
    DEBUG_ASSERT_LE(onv.nbit(), inds.capacity(), "vacant updater inds not large enough for ONV");
    inds.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}


void NdOccupiedUpdater::operator()(const fields::FrmOnv &onv, const defs::inds& map, std::vector<defs::inds> &inds) {
    for (auto& v :inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_dataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            inds[map[ibit]].push_back(ibit);
        }
    }
}

void NdVacantUpdater::operator()(const fields::FrmOnv &onv, const defs::inds& map, std::vector<defs::inds> &inds){
    for (auto& v :inds) v.clear();
    for (size_t idataword = 0ul; idataword < onv.m_dsize; ++idataword) {
        auto work = onv.get_antidataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            inds[map[ibit]].push_back(ibit);
        }
    }
}