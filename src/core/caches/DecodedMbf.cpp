//
// Created by rja on 24/08/2021.
//

#include "DecodedMbf.h"
#include "src/core/basis/AbelianGroup.h"
#include "src/core/field/FrmOnvField.h"
#include "src/core/field/BosOnvField.h"
#include "src/core/field/FrmBosOnvField.h"


void decoded_mbf::SimpleContainer::clear() {
    m_inds.clear();
}

bool decoded_mbf::SimpleContainer::empty() {
    return m_inds.empty();
}

decoded_mbf::frm::Base::Base(const FrmOnvField &mbf) : m_mbf(mbf){}

bool decoded_mbf::frm::Base::is_valid() const {
    return !defs::enable_debug || (m_mbf.hash()==m_last_update_hash);
}

const defs::inds &decoded_mbf::frm::SimpleBase::validated() const {
    DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
    return m_inds;
}

decoded_mbf::frm::SimpleBase::SimpleBase(const FrmOnvField &mbf) : Base(mbf) {}

decoded_mbf::frm::SimpleOccs::SimpleOccs(const FrmOnvField &mbf) : SimpleBase(mbf) {}

const defs::inds& decoded_mbf::frm::SimpleOccs::get() {
    if (!empty()) return validated();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_dataword(idataword);
        while (work) m_inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

decoded_mbf::frm::SimpleVacs::SimpleVacs(const FrmOnvField &mbf) : SimpleBase(mbf) {}

const defs::inds& decoded_mbf::frm::SimpleVacs::get() {
    if (!empty()) return validated();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_antidataword(idataword);
        while (work) m_inds.push_back(bit_utils::next_setbit(work) + idataword * defs::nbit_word);
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

decoded_mbf::frm::LabelledBase::LabelledBase(size_t nelement, const defs::inds &map, const FrmOnvField &mbf) :
        Base(mbf), m_inds(nelement), m_map(map) {
    REQUIRE_LT(*std::max_element(map.cbegin(), map.cend()), nelement,
               "not allocating enough elements in ragged array to accommodate label map");
}

const std::vector<defs::inds> &decoded_mbf::frm::LabelledBase::validated() const {
    DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
    return m_inds;
}

defs::inds decoded_mbf::frm::LabelledBase::make_spinorb_map(const defs::inds &site_irreps, size_t nirrep) {
    auto nsite = site_irreps.size();
    defs::inds out(2*nsite, 0);
    std::copy(site_irreps.cbegin(), site_irreps.cend(), out.begin());
    std::copy(site_irreps.cbegin(), site_irreps.cend(), out.begin()+nsite);
    for (size_t i=nsite; i<nsite*2; ++i) out[i]+=nirrep;
    return out;
}


void decoded_mbf::frm::LabelledBase::clear() {
    for (auto& v: m_inds) v.clear();
    m_simple_inds.clear();
}

bool decoded_mbf::frm::LabelledBase::empty() {
    return m_simple_inds.empty();
}

decoded_mbf::frm::LabelledOccs::LabelledOccs(size_t nelement, const defs::inds &map, const FrmOnvField& mbf) :
        LabelledBase(nelement, map, mbf){}

const std::vector<defs::inds>& decoded_mbf::frm::LabelledOccs::get() {
    if (!empty()) return validated();
    // simple inds are all ready cleared
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_dataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple_inds.push_back(ibit);
        }
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

const defs::inds &decoded_mbf::frm::LabelledOccs::simple() {
    get();
    return m_simple_inds;
}

decoded_mbf::frm::LabelledVacs::LabelledVacs(size_t nelement, const defs::inds &map, const FrmOnvField& mbf) :
        LabelledBase(nelement, map, mbf){}

const std::vector<defs::inds>& decoded_mbf::frm::LabelledVacs::get() {
    if (!empty()) return validated();
    // simple inds are all ready cleared
    for (auto& v : m_inds) v.clear();
    for (size_t idataword = 0ul; idataword < m_mbf.m_dsize; ++idataword) {
        auto work = m_mbf.get_antidataword(idataword);
        while (work) {
            auto ibit = bit_utils::next_setbit(work) + idataword * defs::nbit_word;
            m_inds[m_map[ibit]].push_back(ibit);
            m_simple_inds.push_back(ibit);
        }
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

const defs::inds &decoded_mbf::frm::LabelledVacs::simple() {
    get();
    return m_simple_inds;
}

decoded_mbf::frm::SpinSymOccs::SpinSymOccs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf) :
        NdLabelledOccs<2>({2, grp_map.m_grp.nirrep()},
                          make_spinorb_map(grp_map.m_site_irreps, grp_map.m_grp.nirrep()), mbf) {
}

decoded_mbf::frm::SpinSymVacs::SpinSymVacs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf) :
        NdLabelledVacs<2>({2, grp_map.m_grp.nirrep()},
                          make_spinorb_map(grp_map.m_site_irreps, grp_map.m_grp.nirrep()), mbf) {}

decoded_mbf::frm::NonEmptyPairLabels::NonEmptyPairLabels(
        const FrmOnvField &mbf, decoded_mbf::frm::SpinSymOccs &occs, decoded_mbf::frm::SpinSymVacs &vacs) :
        SimpleBase(mbf), m_occs(occs), m_vacs(vacs) {}

const defs::inds &decoded_mbf::frm::NonEmptyPairLabels::get() {
    if (!empty()) return validated();
    auto &occ = m_occs.get();
    auto &vac = m_vacs.get();
    for (size_t i = 0ul; i < occ.m_format.m_nelement; ++i) {
        if (!occ[i].empty() && !vac[i].empty()) m_inds.push_back(i);
    }
    return m_inds;
}


decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf, const AbelianGroupMap& grp_map):
        m_simple_occs(mbf), m_simple_vacs(mbf),
        m_spin_sym_occs(grp_map, mbf), m_spin_sym_vacs(grp_map, mbf),
        m_nonempty_pair_labels(mbf, m_spin_sym_occs, m_spin_sym_vacs){
    clear();
}

decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf): FrmOnv(mbf, {mbf.m_nsite}){}

void decoded_mbf::FrmOnv::clear() {
    m_simple_occs.clear();
    m_simple_vacs.clear();
    m_spin_sym_occs.clear();
    m_spin_sym_vacs.clear();
    m_nonempty_pair_labels.clear();
    m_occ_sites.clear();
    m_doubly_occ_sites.clear();
    m_not_singly_occ_sites.clear();
}

const defs::inds &decoded_mbf::FrmOnv::occ_sites() {
//    if (m_occ_sites.empty()) {
//        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
//            if (m_mbf.site_nocc(isite)) m_occ_sites.push_back(isite);
//        }
//    }
    return m_occ_sites;
}

const defs::inds &decoded_mbf::FrmOnv::doubly_occ_sites() {
//    if (m_doubly_occ_sites.empty()) {
//        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
//            if (m_mbf.site_nocc(isite)==2) m_doubly_occ_sites.push_back(isite);
//        }
//    }
    return m_doubly_occ_sites;
}

const defs::inds &decoded_mbf::FrmOnv::not_singly_occ_sites() {
//    if (m_not_singly_occ_sites.empty()) {
//        for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
//            if (m_mbf.site_nocc(isite)!=1) m_not_singly_occ_sites.push_back(isite);
//        }
//    }
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
