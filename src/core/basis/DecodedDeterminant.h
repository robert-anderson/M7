//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/FermionOnvSpecifier.h"
#include "src/core/field/Views.h"

/*
 * void updater_fn (const views::FermionOnv&, defs::inds&)
 */

template<typename updater_fn>
class DecodedDeterminant {
    // spin orbital indices
    defs::inds m_inds;

public:
    explicit DecodedDeterminant(const FermionOnvSpecifier &spec) {
        m_inds.reserve(spec.m_nbit);
    }

    explicit DecodedDeterminant(const views::FermionOnv &view) :
            DecodedDeterminant(view.spec()) {
        update(view);
    }

    explicit DecodedDeterminant(const views::FermiBosOnv &view) :
            DecodedDeterminant(view.m_fonv.spec()) {}

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

    void update(const views::FermionOnv &onv) {
        updater_fn()(onv, m_inds);
    };

    void update(const views::FermiBosOnv &onv) {
        updater_fn()(onv.m_fonv, m_inds);
    }
};


struct OccupiedUpdater {
    void operator()(const views::FermionOnv &view, defs::inds &inds);
};


struct VacantUpdater {
    void operator()(const views::FermionOnv &view, defs::inds &inds);
};

typedef DecodedDeterminant<OccupiedUpdater> OccupiedOrbitals;
typedef DecodedDeterminant<VacantUpdater> VacantOrbitals;

#endif //M7_DECODEDDETERMINANT_H
