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

const SpinSymOccOrbs &decoded_mbf::FrmOnv::occ() {
    if (m_occ.empty()) m_occ.update(m_mbf);
    return m_occ;
}

const SpinSymVacOrbs &decoded_mbf::FrmOnv::vac() {
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