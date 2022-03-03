//
// Created by rja on 24/08/2021.
//

#include "CachedOrbs.h"
#include "src/core/basis/AbelianGroup.h"
#include "src/core/field/FrmOnvField.h"
#include "src/core/field/BosOnvField.h"
#include "src/core/field/FrmBosOnvField.h"

decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf, const AbelianGroupMap& grp_map):
        m_mbf(mbf), m_simple_occ(mbf), m_simple_vac(mbf),
        m_spin_sym_occ(grp_map, mbf), m_spin_sym_vac(grp_map, mbf){
    clear();
}

decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf): FrmOnv(mbf, {mbf.m_nsite}){}

void decoded_mbf::FrmOnv::clear() {
    m_simple_occ.clear();
    m_simple_vac.clear();
    m_spin_sym_occ.clear();
    m_spin_sym_vac.clear();
    m_nonempty_pair_labels.clear();
    m_occ_sites.clear();
    m_doubly_occ_sites.clear();
    m_not_singly_occ_sites.clear();
}

const decoded_mbf::spinorbs::SimpleOccs& decoded_mbf::FrmOnv::simple_occs() {
    if (m_simple_occ.empty()) m_simple_occ.update();
    return m_simple_occ;
}

const decoded_mbf::spinorbs::SimpleVacs& decoded_mbf::FrmOnv::simple_vacs() {
    if (m_simple_vac.empty()) m_simple_vac.update();
    return m_simple_vac;
}

const decoded_mbf::spinorbs::SpinSymOccs& decoded_mbf::FrmOnv::spin_sym_occs() {
    if (m_spin_sym_occ.empty()) m_spin_sym_occ.update();
    return m_spin_sym_occ;
}

const decoded_mbf::spinorbs::SpinSymVacs& decoded_mbf::FrmOnv::spin_sym_vacs() {
    if (m_spin_sym_vac.empty()) m_spin_sym_vac.update();
    return m_spin_sym_vac;
}

const defs::inds &decoded_mbf::FrmOnv::nonempty_pair_labels() {
    if (m_nonempty_pair_labels.empty()){
        auto& occ = this->spin_sym_occs();
        auto& vac = this->spin_sym_vacs();
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

void decoded_mbf::spinorbs::SimpleOccs::update() {
    m_inds.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_dataword(idataword);
        while (work) m_inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
}

void decoded_mbf::spinorbs::SimpleVacs::update() {
    m_inds.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_antidataword(idataword);
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

size_t decoded_mbf::spinorbs::LabelledBase::size(const size_t &ielement) const {
    return m_inds[ielement].size();
}

const defs::inds &decoded_mbf::spinorbs::LabelledBase::operator[](const size_t &i) const {
    ASSERT(i<m_inds.size());
    return m_inds[i];
}

void decoded_mbf::spinorbs::LabelledBase::clear() {
    for (auto& v: m_inds) v.clear();
    m_simple_inds.clear();
}

bool decoded_mbf::spinorbs::LabelledBase::empty() {
    return m_simple_inds.empty();
}

decoded_mbf::spinorbs::LabelledOccs::LabelledOccs(size_t nelement, const defs::inds &map, const FrmOnvField& mbf) :
        LabelledBase(nelement, map, mbf){}

void decoded_mbf::spinorbs::LabelledOccs::update() {
    m_simple_inds.clear();
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_dataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple_inds.push_back(ibit);
        }
    }
}

decoded_mbf::spinorbs::LabelledVacs::LabelledVacs(size_t nelement, const defs::inds &map, const FrmOnvField& mbf) :
        LabelledBase(nelement, map, mbf){}

void decoded_mbf::spinorbs::LabelledVacs::update() {
    m_simple_inds.clear();
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_antidataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple_inds.push_back(ibit);
        }
    }
}