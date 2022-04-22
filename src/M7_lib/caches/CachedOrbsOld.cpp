//
// Created by rja on 02/03/2022.
//

#include "CachedOrbsOld.h"


#if 0

CachedOrbs::CachedOrbs(const AbelianGroupMap &grp_map) :
        m_occ(grp_map), m_vac(grp_map){
    clear();
}

CachedOrbs::CachedOrbs(size_t nsite): CachedOrbs(AbelianGroupMap(nsite)){}

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
        for (size_t isite = 0ul; isite < mbf.m_bd.m_nsite; ++isite) {
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
        for (size_t isite = 0ul; isite < mbf.m_bd.m_nsite; ++isite) {
            if (mbf.site_nocc(isite)==2) m_doubly_occ_sites.push_back(isite);
        }
    }
    return m_doubly_occ_sites;
}

const defs::inds &CachedOrbs::not_singly_occ_sites(const field::FrmOnv &mbf) {
    if (m_not_singly_occ_sites.empty()) {
        for (size_t isite = 0ul; isite < mbf.m_bd.m_nsite; ++isite) {
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
#endif