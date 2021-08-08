//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANT_H
#define M7_DECODEDDETERMINANT_H

#include "src/core/field/Fields.h"
#include "AbelianGroup.h"

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

    explicit DecodedDeterminant(const fields::FrmOnv &onv) :
            DecodedDeterminant(onv.nsite()) {
        update(onv);
    }

    explicit DecodedDeterminant(const fields::FrmBosOnv &onv) :
            DecodedDeterminant(onv.nsite()) {
        update(onv);
    }

    DecodedDeterminant(const DecodedDeterminant& other): DecodedDeterminant(other.m_inds.capacity()/2){}
    DecodedDeterminant(DecodedDeterminant&& other): DecodedDeterminant(other.m_inds.capacity()/2){}
    DecodedDeterminant& operator=(const DecodedDeterminant& other){
        m_inds = other.m_inds;
        return *this;
    }
    DecodedDeterminant& operator=(DecodedDeterminant&& other){
        m_inds = std::move(other.m_inds);
        return *this;
    }

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

    void update(const fields::FrmOnv &onv) {
        updater_fn()(onv, m_inds);
    };

    void update(const fields::FrmBosOnv &onv) {
        updater_fn()(onv.m_frm, m_inds);
    }
};


struct OccupiedUpdater {
    void operator()(const fields::FrmOnv &onv, defs::inds &inds);
};


struct VacantUpdater {
    void operator()(const fields::FrmOnv &onv, defs::inds &inds);
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
        m_format(shape), m_inds(m_format.m_nelement), m_map(map){
        for (auto& v: m_inds) v.reserve(2*nsite);
        ASSERT(m_map.size()==2*nsite);
        ASSERT(*std::max_element(map.cbegin(), map.cend())<m_format.m_nelement);
    }

    NdDecodedDeterminant(std::array<size_t, nind> shape, const fields::FrmOnv &onv, const defs::inds& map) :
            NdDecodedDeterminant(shape, onv.m_nsite, map) {
        update(onv);
    }

    NdDecodedDeterminant(std::array<size_t, nind> shape, const fields::FrmBosOnv &onv, const defs::inds& map) :
            NdDecodedDeterminant(shape, onv.m_frm.m_nsite, map) {}


    NdDecodedDeterminant(const NdDecodedDeterminant& other):
        NdDecodedDeterminant(other.m_format.m_shape, other.m_inds.capacity()/2, other.m_map){}
    NdDecodedDeterminant(NdDecodedDeterminant&& other):
            NdDecodedDeterminant(other.m_format.m_shape, other.m_inds.capacity()/2, other.m_map){}
    NdDecodedDeterminant& operator=(const NdDecodedDeterminant& other){
        m_inds = other.m_inds;
        return *this;
    }
    NdDecodedDeterminant& operator=(NdDecodedDeterminant&& other){
        m_inds = std::move(other.m_inds);
        return *this;
    }

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

    void update(const fields::FrmOnv &onv) {
        updater_fn()(onv, m_map, m_inds);
    };

    void update(const fields::FrmBosOnv &onv) {
        updater_fn()(onv.m_frm, m_map, m_inds);
    }
};


struct NdOccupiedUpdater {
    void operator()(const fields::FrmOnv &onv, const defs::inds& map, std::vector<defs::inds> &inds);
};


struct NdVacantUpdater {
    void operator()(const fields::FrmOnv &onv, const defs::inds& map, std::vector<defs::inds> &inds);
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

template<typename updater_fn>
class SymmDecodedDeterminant : public NdDecodedDeterminant<updater_fn, 2> {

    defs::inds make_map(const AbelianGroupMap& grp_map) {
        defs::inds out(2*grp_map.m_nsite, 0);
        std::copy(grp_map.m_site_irreps.cbegin(), grp_map.m_site_irreps.cend(), out.begin());
        std::copy(grp_map.m_site_irreps.cbegin(), grp_map.m_site_irreps.cend(), out.begin()+grp_map.m_nsite);
        for (size_t i=grp_map.m_nsite; i<grp_map.m_nsite*2; ++i) out[i]+=grp_map.m_grp.nirrep();
        return out;
    }

public:
    explicit SymmDecodedDeterminant(const AbelianGroupMap& grp_map) :
            NdDecodedDeterminant<updater_fn, 2>({2, grp_map.m_grp.nirrep()}, grp_map.m_nsite, make_map(grp_map)){
    }
};

typedef SymmDecodedDeterminant<NdOccupiedUpdater> SymmOccupiedOrbitals;
typedef SymmDecodedDeterminant<NdVacantUpdater> SymmVacantOrbitals;


#endif //M7_DECODEDDETERMINANT_H
