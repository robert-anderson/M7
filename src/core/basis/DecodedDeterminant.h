//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/Fields.h"

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

    explicit DecodedDeterminant(const fields::Onv<0> &onv) :
            DecodedDeterminant(onv.m_nsite) {
        update(onv);
    }

    explicit DecodedDeterminant(const fields::Onv<1> &onv) :
            DecodedDeterminant(onv.m_frm.m_nsite) {}

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

    void update(const fields::Onv<0> &onv) {
        updater_fn()(onv, m_inds);
    };

    void update(const fields::Onv<1> &onv) {
        updater_fn()(onv.m_frm, m_inds);
    }
};


struct OccupiedUpdater {
    void operator()(const fields::Onv<0> &view, defs::inds &inds);
};


struct VacantUpdater {
    void operator()(const fields::Onv<0> &view, defs::inds &inds);
};

typedef DecodedDeterminant<OccupiedUpdater> OccupiedOrbitals;
typedef DecodedDeterminant<VacantUpdater> VacantOrbitals;



template<typename updater_fn, size_t nind>
class NdDecodedDeterminant {
    const NdFormat<nind> m_format;
    std::vector<defs::inds> m_inds;
    const defs::inds m_map;

public:
    NdDecodedDeterminant(std::array<size_t, nind> shape, size_t nsite, const defs::inds& map):
        m_format(shape), m_inds(m_format.nelement()), m_map(map){
        for (auto& v: m_inds) v.reserve(2*nsite);
        ASSERT(m_map.size()==2*nsite);
        ASSERT(*std::max_element(map.cbegin(), map.cend())<m_format.nelement());
    }

    NdDecodedDeterminant(std::array<size_t, nind> shape, const fields::Onv<0> &onv, const defs::inds& map) :
            NdDecodedDeterminant(shape, onv.m_nsite, map) {
        update(onv);
    }

    NdDecodedDeterminant(std::array<size_t, nind> shape, const fields::Onv<1> &onv, const defs::inds& map) :
            NdDecodedDeterminant(shape, onv.m_frm.m_nsite, map) {}

    size_t size(const size_t& ielement) const {
        return m_inds[ielement].size();
    }

    size_t size(const std::array<size_t, nind>& inds) const {
        return m_inds[m_format.flatten(inds)].size();
    }

    const defs::inds& operator[](const size_t& i) const{
        ASSERT(i<m_inds.size());
        return m_inds[i];
    }

    const defs::inds& operator[](const std::array<size_t, nind>& inds) const{
        return m_inds[m_format.flatten(inds)];
    }

    void update(const fields::Onv<0> &onv) {
        updater_fn()(onv, m_map, m_inds);
    };

    void update(const fields::Onv<1> &onv) {
        updater_fn()(onv.m_frm, m_map, m_inds);
    }
};


struct NdOccupiedUpdater {
    void operator()(const fields::Onv<0> &view, const defs::inds& map, std::vector<defs::inds> &inds);
};


struct NdVacantUpdater {
    void operator()(const fields::Onv<0> &view, const defs::inds& map, std::vector<defs::inds> &inds);
};

template<size_t nind>
using NdOccupiedOrbitals = NdDecodedDeterminant<NdOccupiedUpdater, nind>;
template<size_t nind>
using NdVacantOrbitals = NdDecodedDeterminant<NdVacantUpdater, nind>;

template<typename updater_fn>
class SpinDecodedDeterminant : public NdDecodedDeterminant<updater_fn, 1> {

    defs::inds make_spin_map(size_t nsite) {
        defs::inds out(2*nsite, 0);
        for (auto it = out.begin()+nsite; it!=out.end(); ++it) *it = 1ul;
        return out;
    }

public:
    explicit SpinDecodedDeterminant(size_t nsite) :
        NdDecodedDeterminant<updater_fn, 1>({2}, nsite, make_spin_map(nsite)){}
};

typedef SpinDecodedDeterminant<NdOccupiedUpdater> SpinOccupiedOrbitals;
typedef SpinDecodedDeterminant<NdVacantUpdater> SpinVacantOrbitals;

#endif //M7_DECODEDDETERMINANT_H
