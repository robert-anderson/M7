//
// Created by rja on 24/08/2021.
//

#include "CachedOrbs.h"
#include "src/core/field/FrmOnvField.h"
#include "src/core/field/BosOnvField.h"
#include "src/core/field/FrmBosOnvField.h"

decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf, const AbelianGroupMap& grp_map):
        m_mbf(mbf), m_occ(grp_map), m_vac(grp_map){
    clear();
}

void decoded_mbf::FrmOnv::clear() {
    m_occ.clear();
    m_vac.clear();
    m_nonempty_pair_labels.clear();
    m_occ_sites.clear();
    m_doubly_occ_sites.clear();
    m_not_singly_occ_sites.clear();
}

const decoded_mbf::spinorbs::SpinSymOccs& decoded_mbf::FrmOnv::occ() {
    if (m_occ.empty()) m_occ.update(m_mbf);
    return m_occ;
}

const decoded_mbf::spinorbs::SpinSymVacs& decoded_mbf::FrmOnv::vac() {
    if (m_vac.empty()) m_vac.update(m_mbf);
    return m_vac;
}

const defs::inds &decoded_mbf::FrmOnv::nonempty_pair_labels() {
    if (m_nonempty_pair_labels.empty()){
        auto& occ = this->occ();
        auto& vac = this->vac();
        for (size_t i=0ul; i<occ.m_format.m_nelement; ++i){
            if (!occ[i].empty() && !vac[i].empty()) m_nonempty_pair_labels.push_back(i);
        }
    }
    return m_nonempty_pair_labels;
}

const defs::inds &decoded_mbf::FrmOnv::occ_sites() {
    if (m_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
            if (m_mbf.site_nocc(isite)) m_occ_sites.push_back(isite);
        }
    }
    return m_occ_sites;
}

const defs::inds &decoded_mbf::FrmOnv::doubly_occ_sites() {
    if (m_doubly_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
            if (m_mbf.site_nocc(isite)==2) m_doubly_occ_sites.push_back(isite);
        }
    }
    return m_doubly_occ_sites;
}

const defs::inds &decoded_mbf::FrmOnv::not_singly_occ_sites() {
    if (m_not_singly_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
            if (m_mbf.site_nocc(isite)!=1) m_not_singly_occ_sites.push_back(isite);
        }
    }
    return m_not_singly_occ_sites;
}


decoded_mbf::BosOnv::BosOnv(const BosOnvField &mbf) : m_mbf(mbf){}

void decoded_mbf::BosOnv::clear() {
    m_bos_op_inds.clear();
    m_occ_bos_inds.clear();
}

const defs::inds &decoded_mbf::BosOnv::bos_op_inds() {
    if (m_bos_op_inds.empty()) {
        for (size_t imode = 0ul; imode < m_mbf.m_nelement; ++imode) {
            for (size_t iop =0ul; iop<m_mbf[imode]; ++iop){
                m_bos_op_inds.push_back(imode);
            }
        }
    }
    return m_bos_op_inds;
}

const defs::inds &decoded_mbf::BosOnv::occ_bos_inds() {
    if (m_occ_bos_inds.empty()) {
        for (size_t imode = 0ul; imode < m_mbf.m_nelement; ++imode) {
            if (m_mbf[imode]) m_occ_bos_inds.push_back(imode);
        }
    }
    return m_occ_bos_inds;
}

decoded_mbf::FrmBosOnv::FrmBosOnv(const FrmBosOnvField &mbf) : m_mbf(mbf){}

void decoded_mbf::FrmBosOnv::clear() {
    m_occ_sites_nonzero_bosons.clear(); // frm-bos cached assets must cleared in both frm and bos clear methods
}

const defs::inds &decoded_mbf::FrmBosOnv::occ_sites_nonzero_bosons() {
    if (m_occ_sites_nonzero_bosons.empty()) {
        // TODO
//        const auto& occ = occ_sites(m_mbf.m_frm);
//        for (auto& imode: occ) if (mbf.m_bos[imode]) m_occ_sites_nonzero_bosons.push_back(imode);
    }
    return m_occ_sites_nonzero_bosons;
}



CachedOrbs::CachedOrbs(const AbelianGroupMap &grp_map) :
        m_occ(grp_map), m_vac(grp_map){
    clear();
}

void CachedOrbs::clear_frmbos_only() {
    m_occ_sites_nonzero_bosons.clear(); // frm-bos cached assets must cleared in both frm and bos clear methods
}

void CachedOrbs::clear_frm() {
    m_occ.clear();
    m_vac.clear();
    m_nonempty_pair_labels.clear();
    m_occ_sites.clear();
    m_doubly_occ_sites.clear();
    m_not_singly_occ_sites.clear();
    clear_frmbos_only();
}

void CachedOrbs::clear(const field::FrmOnv& mbf) {
    clear_frm();
}

void CachedOrbs::clear_bos() {
    m_bos_op_inds.clear();
    m_occ_bos_inds.clear();
    clear_frmbos_only();
}

void CachedOrbs::clear(const field::BosOnv& mbf) {
    clear_bos();
}

void CachedOrbs::clear(const field::FrmBosOnv& mbf) {
    clear(mbf.m_frm);
    clear(mbf.m_bos);
}

void CachedOrbs::clear() {
    clear_frm();
    clear_bos();
}

const SpinSymOccOrbs__ &CachedOrbs::occ(const field::FrmOnv &mbf) {
    if (m_occ.empty()) m_occ.update(mbf);
    return m_occ;
}

const SpinSymVacOrbs__ &CachedOrbs::vac(const field::FrmOnv &mbf) {
    if (m_vac.empty()) m_vac.update(mbf);
    return m_vac;
}

const defs::inds &CachedOrbs::nonempty_pair_labels(const field::FrmOnv &mbf) {
    if (m_nonempty_pair_labels.empty()){
        auto& occ = this->occ(mbf);
        auto& vac = this->vac(mbf);
        for (size_t i=0ul; i<occ.m_format.m_nelement; ++i){
            if (!occ[i].empty() && !vac[i].empty()) m_nonempty_pair_labels.push_back(i);
        }
    }
    return m_nonempty_pair_labels;
}

const defs::inds &CachedOrbs::occ_sites(const field::FrmOnv &mbf) {
    if (m_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < mbf.nsite(); ++isite) {
            if (mbf.site_nocc(isite)) m_occ_sites.push_back(isite);
        }
    }
    return m_occ_sites;
}

const defs::inds &CachedOrbs::occ_sites_nonzero_bosons(const field::FrmBosOnv &mbf) {
    if (m_occ_sites_nonzero_bosons.empty()) {
        const auto& occ = occ_sites(mbf.m_frm);
        for (auto& imode: occ) if (mbf.m_bos[imode]) m_occ_sites_nonzero_bosons.push_back(imode);
    }
    return m_occ_sites_nonzero_bosons;
}

const defs::inds &CachedOrbs::doubly_occ_sites(const field::FrmOnv &mbf) {
    if (m_doubly_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < mbf.nsite(); ++isite) {
            if (mbf.site_nocc(isite)==2) m_doubly_occ_sites.push_back(isite);
        }
    }
    return m_doubly_occ_sites;
}

const defs::inds &CachedOrbs::not_singly_occ_sites(const field::FrmOnv &mbf) {
    if (m_not_singly_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < mbf.nsite(); ++isite) {
            if (mbf.site_nocc(isite)!=1) m_not_singly_occ_sites.push_back(isite);
        }
    }
    return m_not_singly_occ_sites;
}

const defs::inds &CachedOrbs::bos_op_inds(const field::BosOnv &mbf) {
    if (m_bos_op_inds.empty()) {
        for (size_t imode = 0ul; imode < mbf.m_nelement; ++imode) {
            for (size_t iop =0ul; iop<mbf[imode]; ++iop){
                m_bos_op_inds.push_back(imode);
            }
        }
    }
    return m_bos_op_inds;
}

const defs::inds &CachedOrbs::occ_bos_inds(const field::BosOnv &mbf) {
    if (m_occ_bos_inds.empty()) {
        for (size_t imode = 0ul; imode < mbf.m_nelement; ++imode) {
            if (mbf[imode]) m_occ_bos_inds.push_back(imode);
        }
    }
    return m_occ_bos_inds;
}











size_t decoded_mbf::spinorbs::SimpleBase::size() const {
    return m_inds.size();
}

const size_t &decoded_mbf::spinorbs::SimpleBase::operator[](const size_t &i) const {
    ASSERT(i<size());
    return m_inds[i];
}

const defs::inds &decoded_mbf::spinorbs::SimpleBase::inds() const {
    return m_inds;
}

void decoded_mbf::spinorbs::SimpleBase::clear() {
    m_inds.clear();
}

bool decoded_mbf::spinorbs::SimpleBase::empty() {
    return m_inds.empty();
}

void decoded_mbf::spinorbs::SimpleOccs::update(const FrmOnvField &mbf) {
    m_inds.clear();
    for (size_t idataword = 0ul; idataword < mbf.m_dsize; ++idataword) {
        auto work = mbf.get_dataword(idataword);
        while (work) m_inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void decoded_mbf::spinorbs::SimpleVacs::update(const FrmOnvField &mbf) {
    m_inds.clear();
    for (size_t idataword = 0ul; idataword < mbf.m_dsize; ++idataword) {
        auto work = mbf.get_antidataword(idataword);
        while (work) m_inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

defs::inds decoded_mbf::spinorbs::LabelledBase::make_spinorb_map(const defs::inds &site_irreps, size_t nirrep) {
    auto nsite = site_irreps.size();
    defs::inds out(2*nsite, 0);
    std::copy(site_irreps.cbegin(), site_irreps.cend(), out.begin());
    std::copy(site_irreps.cbegin(), site_irreps.cend(), out.begin()+nsite);
    for (size_t i=nsite; i<nsite*2; ++i) out[i]+=nirrep;
    return out;
}

decoded_mbf::spinorbs::LabelledBase::LabelledBase(size_t nelement, const defs::inds &map) :
    m_inds(nelement), m_map(map){
    REQUIRE_LT(*std::max_element(map.cbegin(), map.cend()), nelement,
               "not allocating enough elements in ragged array to accommodate label map");
}

decoded_mbf::spinorbs::LabelledBase::LabelledBase(const decoded_mbf::spinorbs::LabelledBase &other) :
        LabelledBase(other.m_inds.size(), other.m_map){}

decoded_mbf::spinorbs::LabelledBase::LabelledBase(decoded_mbf::spinorbs::LabelledBase &&other) :
        LabelledBase(other.m_inds.size(), other.m_map){}

decoded_mbf::spinorbs::LabelledBase &decoded_mbf::spinorbs::LabelledBase::operator=(
        const decoded_mbf::spinorbs::LabelledBase &other) {
    DEBUG_ASSERT_EQ(m_map, other.m_map, "label maps do not match");
    m_inds = other.m_inds;
    return *this;
}

decoded_mbf::spinorbs::LabelledBase &decoded_mbf::spinorbs::LabelledBase::operator=(
        decoded_mbf::spinorbs::LabelledBase &&other) {
    DEBUG_ASSERT_EQ(m_map, other.m_map, "label maps do not match");
    m_inds = std::move(other.m_inds);
    return *this;
}

size_t decoded_mbf::spinorbs::LabelledBase::size(const size_t &ielement) const {
    return m_inds[ielement].size();
}

const defs::inds &decoded_mbf::spinorbs::LabelledBase::operator[](const size_t &i) const {
    ASSERT(i<m_inds.size());
    return m_inds[i];
}

void decoded_mbf::spinorbs::LabelledBase::clear() {
    for (auto& v: m_inds) v.clear();
}

decoded_mbf::spinorbs::LabelledOccs::LabelledOccs(size_t nelement, const defs::inds &map) :
        LabelledBase(nelement, map){}

void decoded_mbf::spinorbs::LabelledOccs::update(const FrmOnvField &mbf) {
    m_simple.clear();
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < mbf.m_dsize; ++idataword) {
        auto work = mbf.get_dataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple.m_inds.push_back(ibit);
        }
    }
}

void decoded_mbf::spinorbs::LabelledOccs::clear() {
    LabelledBase::clear();
    m_simple.clear();
}

bool decoded_mbf::spinorbs::LabelledOccs::empty() {
    return m_simple.empty();
}

decoded_mbf::spinorbs::LabelledVacs::LabelledVacs(size_t nelement, const defs::inds &map) :
        LabelledBase(nelement, map){}

void decoded_mbf::spinorbs::LabelledVacs::update(const FrmOnvField &mbf) {
    m_simple.clear();
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < mbf.m_dsize; ++idataword) {
        auto work = mbf.get_antidataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple.m_inds.push_back(ibit);
        }
    }
}

void decoded_mbf::spinorbs::LabelledVacs::clear() {
    LabelledBase::clear();
    m_simple.clear();
}

bool decoded_mbf::spinorbs::LabelledVacs::empty() {
    return m_simple.empty();
}