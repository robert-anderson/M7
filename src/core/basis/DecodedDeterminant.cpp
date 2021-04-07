//
// Created by Robert John Anderson on 2020-03-30.
//

#include "DecodedDeterminant.h"

void OccupiedUpdater::operator()(const fields::Onv<0> &view, defs::inds &inds) {
    ASSERT(view.nbit() <= inds.capacity());
    inds.clear();
    for (size_t idataword = 0ul; idataword < view.m_dsize; ++idataword) {
        auto work = view.get_dataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
    }
}

void VacantUpdater::operator()(const fields::Onv<0> &view, defs::inds &inds) {
    ASSERT(view.nbit() <= inds.capacity());
    inds.clear();
    for (size_t idataword = 0ul; idataword < view.m_dsize; ++idataword) {
        auto work = view.get_antidataword(idataword);
        while (work) inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_data);
    }
}
