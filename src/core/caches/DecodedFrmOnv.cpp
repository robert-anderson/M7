//
// Created by rja on 04/03/2022.
//

#include "DecodedFrmOnv.h"
#include "src/core/field/FrmOnvField.h"

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

decoded_mbf::frm::NonEmptyPairLabels::NonEmptyPairLabels(const FrmOnvField &mbf) : SimpleBase(mbf){}

const defs::inds &decoded_mbf::frm::NonEmptyPairLabels::get() {
    if (!empty()) return validated();
    auto &occ = m_mbf.m_decoded.m_spin_sym_occs.get();
    auto &vac = m_mbf.m_decoded.m_spin_sym_vacs.get();
    for (size_t i = 0ul; i < occ.m_format.m_nelement; ++i) {
        if (!occ[i].empty() && !vac[i].empty()) m_inds.push_back(i);
    }
    return m_inds;
}

decoded_mbf::frm::OccSites::OccSites(const FrmOnvField &mbf) : SimpleBase(mbf){}

const defs::inds &decoded_mbf::frm::OccSites::get() {
    if (!empty()) return validated();
    for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
        if (m_mbf.site_nocc(isite)) m_inds.push_back(isite);
    }
    return m_inds;
}

decoded_mbf::frm::DoublyOccSites::DoublyOccSites(const FrmOnvField &mbf) : SimpleBase(mbf){}

const defs::inds &decoded_mbf::frm::DoublyOccSites::get() {
    if (!empty()) return validated();
    for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
        if (m_mbf.site_nocc(isite)==2) m_inds.push_back(isite);
    }
    return m_inds;
}

decoded_mbf::frm::NotSinglyOccSites::NotSinglyOccSites(const FrmOnvField &mbf) : SimpleBase(mbf){}

const defs::inds &decoded_mbf::frm::NotSinglyOccSites::get() {
    if (!empty()) return validated();
    for (size_t isite = 0ul; isite < m_mbf.nsite(); ++isite) {
        if (m_mbf.site_nocc(isite)!=1) m_inds.push_back(isite);
    }
    return m_inds;
}



decoded_mbf::FrmOnv::FrmOnv(const FrmOnvField& mbf, const AbelianGroupMap& grp_map):
        m_simple_occs(mbf), m_simple_vacs(mbf),
        m_spin_sym_occs(grp_map, mbf), m_spin_sym_vacs(grp_map, mbf),
        m_nonempty_pair_labels(mbf), m_occ_sites(mbf),
        m_doubly_occ_sites(mbf), m_not_singly_occ_sites(mbf){
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