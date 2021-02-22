//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/fieldz/FieldsZ.h"

/*
 * void updater_fn (const views::FermionOnv&, defs::inds&)
 */

template<typename updater_fn>
class DecodedDeterminant {
    // spin orbital indices
    defs::inds m_inds;

public:
    explicit DecodedDeterminant(size_t nsite) {
        m_inds.reserve(2*nsite);
    }

    explicit DecodedDeterminant(const fieldsz::Onv<0> &onv) :
            DecodedDeterminant(onv.m_nsite) {
        update(onv);
    }

    explicit DecodedDeterminant(const fieldsz::Onv<1> &onv) :
            DecodedDeterminant(onv.m_fonv.m_nsite) {}

    size_t size() const {
        return m_inds.size();
    }

    const size_t& operator[](const size_t& i) const{
        ASSERT(i<size());
        return m_inds[i];
    }

    const defs::inds& inds() const{
        return m_inds;
    }

    void update(const fieldsz::Onv<0> &onv) {
        updater_fn()(onv, m_inds);
    };

    void update(const fieldsz::Onv<1> &onv) {
        updater_fn()(onv.m_fonv, m_inds);
    }
};


struct OccupiedUpdater {
    void operator()(const fieldsz::Onv<0> &view, defs::inds &inds);
};


struct VacantUpdater {
    void operator()(const fieldsz::Onv<0> &view, defs::inds &inds);
};

typedef DecodedDeterminant<OccupiedUpdater> OccupiedOrbitals;
typedef DecodedDeterminant<VacantUpdater> VacantOrbitals;

#endif //M7_DECODEDDETERMINANT_H
