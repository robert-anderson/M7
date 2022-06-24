//
// Created by Robert J. Anderson on 04/03/2022.
//

#include <M7_lib/field/FrmBosOnvField.h>

#include "DecodedFrmBosOnv.h"

decoded_mbf::frmbos::Base::Base(const FrmBosOnvField &mbf) : m_mbf(mbf){}

bool decoded_mbf::frmbos::Base::is_valid() const {
    return !c_enable_debug || (m_mbf.hash() == m_last_update_hash);
}

const uintv_t &decoded_mbf::frmbos::SimpleBase::validated() const {
    DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
    return m_inds;
}

decoded_mbf::frmbos::SimpleBase::SimpleBase(const FrmBosOnvField &mbf) : Base(mbf){}

decoded_mbf::frmbos::OccSitesNonzeroBosons::OccSitesNonzeroBosons(const FrmBosOnvField &mbf) : SimpleBase(mbf){}

const uintv_t &decoded_mbf::frmbos::OccSitesNonzeroBosons::get() {
    DEBUG_ASSERT_EQ(m_mbf.m_frm.m_basis.m_nsite, m_mbf.m_bos.m_basis.m_nmode,
                    "this cache only makes sense with a 1-to-1 correspondence between modes and sites");
    if (!empty()) return validated();
    const auto& occ = m_mbf.m_frm.m_decoded.m_occ_sites.get();
    for (auto& imode: occ) if (m_mbf.m_bos[imode]) m_inds.push_back(imode);
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

decoded_mbf::FrmBosOnv::FrmBosOnv(const FrmBosOnvField &mbf) : m_mbf(mbf), m_occ_sites_nonzero_bosons(mbf){}

void decoded_mbf::FrmBosOnv::clear() {
    m_occ_sites_nonzero_bosons.clear();
    // frm-bos cached assets must cleared in both frm and bos clear methods
    m_mbf.m_frm.m_decoded.clear();
    m_mbf.m_bos.m_decoded.clear();
}
