//
// Created by Robert John Anderson on 2020-03-30.
//

#ifndef M7_DECODEDDETERMINANTS_H
#define M7_DECODEDDETERMINANTS_H

#include "src/core/field/Fields.h"
#include "AbelianGroup.h"

/*
 * void updater_fn (const views::FermionOnv&, defs::inds&)
 */

template<typename updater_fn>
struct FlatOrbs {
    /**
     * spin orbital indices
     */
    defs::inds m_inds;

public:
    explicit FlatOrbs(size_t nsite) {
        m_inds.reserve(2*nsite);
    }

    explicit FlatOrbs(const field::FrmOnv &mbf) :
            FlatOrbs(mbf.nsite()) {
        update(mbf);
    }

    explicit FlatOrbs(const field::FrmBosOnv &mbf) :
            FlatOrbs(mbf.nsite()) {
        update(mbf);
    }

    FlatOrbs(const FlatOrbs& other): FlatOrbs(other.m_inds.capacity() / 2){}
    FlatOrbs(FlatOrbs&& other): FlatOrbs(other.m_inds.capacity() / 2){}
    FlatOrbs& operator=(const FlatOrbs& other){
        m_inds = other.m_inds;
        return *this;
    }
    FlatOrbs& operator=(FlatOrbs&& other){
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

    void update(const field::FrmOnv &mbf) {
        updater_fn()(mbf, m_inds);
    };

    void update(const field::FrmBosOnv &mbf) {
        updater_fn()(mbf.m_frm, m_inds);
    }

    void clear() {
        m_inds.clear();
    }
    bool empty() {
        return m_inds.empty();
    }
};


struct OccupiedUpdater {
    void operator()(const field::FrmOnv &mbf, defs::inds &inds);
};


struct VacantUpdater {
    void operator()(const field::FrmOnv &mbf, defs::inds &inds);
};

typedef FlatOrbs<OccupiedUpdater> OccOrbs;
typedef FlatOrbs<VacantUpdater> VacOrbs;


template<typename updater_fn, size_t nind>
struct NdOrbs {
    const NdFormat<nind> m_format;
    std::vector<defs::inds> m_inds;
    const defs::inds m_map;

    FlatOrbs<updater_fn> m_flat;
    NdOrbs(std::array<size_t, nind> shape, size_t nsite, const defs::inds& map):
        m_format(shape), m_inds(m_format.m_nelement), m_map(map), m_flat(nsite){
        if (!nsite) return;
        for (auto& v: m_inds) v.reserve(2*nsite);
        ASSERT(m_map.size()==2*nsite);
        ASSERT(*std::max_element(map.cbegin(), map.cend())<m_format.m_nelement);
    }

    NdOrbs(std::array<size_t, nind> shape, const field::FrmOnv &mbf, const defs::inds& map) :
            NdOrbs(shape, mbf.m_nsite, map) {
        update(mbf);
    }

    NdOrbs(std::array<size_t, nind> shape, const field::FrmBosOnv &mbf, const defs::inds& map) :
            NdOrbs(shape, mbf.m_frm.m_nsite, map) {}


    NdOrbs(const NdOrbs& other):
            NdOrbs(other.m_format.m_shape, other.m_inds.capacity() / 2, other.m_map){}
    NdOrbs(NdOrbs&& other):
            NdOrbs(other.m_format.m_shape, other.m_inds.capacity() / 2, other.m_map){}
    NdOrbs& operator=(const NdOrbs& other){
        m_inds = other.m_inds;
        return *this;
    }
    NdOrbs& operator=(NdOrbs&& other){
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

    void update(const field::FrmOnv &mbf) {
        updater_fn()(mbf, m_map, m_flat.m_inds, m_inds);
    };

    void update(const field::FrmBosOnv &mbf) {
        updater_fn()(mbf.m_frm, m_map, m_flat.m_inds, m_inds);
    }

    /**
     * doesn't actually clear the Nd partitions, just the flat. the update calls always clear m_inds anyway
     */
    void clear() {
        m_flat.clear();
    }
    bool empty() {
        return m_flat.empty();
    }
};


struct NdOccupiedUpdater {
    /**
     * update the flat and N-dimensional ragged arrays of set bits in the fermion ONV
     * @param mbf
     *  many fermion basis function to decode
     * @param map
     *  map between spin orbitals and their label (i.e. element of nd_inds they append to)
     * @param flat_inds
     *  the vector storing set bits of any label
     * @param nd_inds
     *  the vector of vectors storing set bits in each label
     */
    void operator()(const field::FrmOnv &mbf, const defs::inds& map,
            defs::inds& flat_inds, std::vector<defs::inds> &nd_inds);
};


struct NdVacantUpdater {
    /**
     * update the flat and N-dimensional ragged arrays of clear bits in the fermion ONV
     * @param mbf
     *  many fermion basis function to decode
     * @param map
     *  map between spin orbitals and their label (i.e. element of nd_inds they append to)
     * @param flat_inds
     *  the vector storing clear bits of any label
     * @param nd_inds
     *  the vector of vectors storing clear bits in each label
     */
    void operator()(const field::FrmOnv &mbf, const defs::inds& map,
            defs::inds& flat_inds, std::vector<defs::inds> &nd_inds);
};

template<size_t nind>
using NdOccOrbs = NdOrbs<NdOccupiedUpdater, nind>;
template<size_t nind>
using NdVacOrbs = NdOrbs<NdVacantUpdater, nind>;

template<typename updater_fn>
class SpinNdOrbs : public NdOrbs<updater_fn, 1> {
    defs::inds make_spin_map(size_t nsite) {
        defs::inds out(2*nsite, 0);
        for (auto it = out.begin()+nsite; it!=out.end(); ++it) *it = 1ul;
        return out;
    }

public:
    explicit SpinNdOrbs(size_t nsite) :
            NdOrbs<updater_fn, 1>({2}, nsite, make_spin_map(nsite)){}
};

typedef SpinNdOrbs<NdOccupiedUpdater> SpinOccOrbs;
typedef SpinNdOrbs<NdVacantUpdater> SpinVacOrbs;

template<typename updater_fn>
class SpinSymNdOrbs : public NdOrbs<updater_fn, 2> {

    defs::inds make_map(const AbelianGroupMap& grp_map) {
        defs::inds out(2*grp_map.m_nsite, 0);
        std::copy(grp_map.m_site_irreps.cbegin(), grp_map.m_site_irreps.cend(), out.begin());
        std::copy(grp_map.m_site_irreps.cbegin(), grp_map.m_site_irreps.cend(), out.begin()+grp_map.m_nsite);
        for (size_t i=grp_map.m_nsite; i<grp_map.m_nsite*2; ++i) out[i]+=grp_map.m_grp.nirrep();
        return out;
    }

public:
    explicit SpinSymNdOrbs(const AbelianGroupMap& grp_map) :
            NdOrbs<updater_fn, 2>({2, grp_map.m_grp.nirrep()}, grp_map.m_nsite, make_map(grp_map)){
    }
};

typedef SpinSymNdOrbs<NdOccupiedUpdater> SpinSymOccOrbs;
typedef SpinSymNdOrbs<NdVacantUpdater> SpinSymVacOrbs;


#endif //M7_DECODEDDETERMINANTS_H
