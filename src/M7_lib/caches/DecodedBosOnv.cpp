//
// Created by Robert J. Anderson on 04/03/2022.
//

#include <M7_lib/field/BosOnvField.h>

#include "DecodedBosOnv.h"

decoded_mbf::bos::Base::Base(const BosOnvField &mbf) : m_mbf(mbf){}

bool decoded_mbf::bos::Base::is_valid() const {
    return !defs::enable_debug || (m_mbf.hash()==m_last_update_hash);
}

const defs::inds_t &decoded_mbf::bos::SimpleBase::validated() const {
    DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
    return m_inds;
}

decoded_mbf::bos::SimpleBase::SimpleBase(const BosOnvField &mbf) : Base(mbf){}

decoded_mbf::bos::Expanded::Expanded(const BosOnvField &mbf) : SimpleBase(mbf){}

const defs::inds_t &decoded_mbf::bos::Expanded::get() {
    if (!empty()) return validated();
    for (size_t imode = 0ul; imode < m_mbf.m_nelement; ++imode) {
        for (size_t iop = 0ul; iop < m_mbf[imode]; ++iop) {
            m_inds.push_back(imode);
        }
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

decoded_mbf::bos::OccModes::OccModes(const BosOnvField &mbf) : SimpleBase(mbf){}

const defs::inds_t &decoded_mbf::bos::OccModes::get() {
    if (!empty()) return validated();
    for (size_t imode = 0ul; imode < m_mbf.m_nelement; ++imode) {
        if (m_mbf[imode]) m_inds.push_back(imode);
    }
    m_last_update_hash = m_mbf.hash();
    return m_inds;
}

decoded_mbf::BosOnv::BosOnv(const BosOnvField &mbf) : m_expanded(mbf), m_occ_modes(mbf){}

void decoded_mbf::BosOnv::clear() {
    m_expanded.clear();
    m_occ_modes.clear();
}
