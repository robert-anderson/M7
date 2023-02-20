//
// Created by Robert J. Anderson on 12/08/2021.
//

#include "RdmIndsField.h"

uinta_t<4> RdmIndsField::make_nops() const {
    return {m_exsig.nfrm_cre(), m_exsig.nfrm_ann(), m_exsig.nbos_cre(), m_exsig.nbos_ann()};
}

uinta_t<4> RdmIndsField::make_nop_offsets() const {
    uinta_t<4> nop_offsets{};
    for (uint_t i = 1ul; i < 4; ++i) nop_offsets[i] = nop_offsets[i - 1] + m_nops[i - 1];
    return nop_offsets;
}

RdmIndsField::RdmIndsField(Row *row, OpSig exsig, str_t name) :
    NdNumberField<mae_ind_t, 1>(row, {exsig.nop()}, name),
    m_exsig(exsig), m_nops(make_nops()), m_nop_offsets(make_nop_offsets()),
    m_frm(*this, m_nop_offsets[0], m_nops[0], m_nop_offsets[1], m_nops[1]),
    m_bos(*this, m_nop_offsets[2], m_nops[2], m_nop_offsets[3], m_nops[3]){}

RdmIndsField::RdmIndsField(const RdmIndsField &other) :
        NdNumberField<mae_ind_t, 1>(other), m_exsig(other.m_exsig), m_nops(make_nops()), m_nop_offsets(make_nop_offsets()),
        m_frm(*this, m_nop_offsets[0], m_nops[0], m_nop_offsets[1], m_nops[1]),
        m_bos(*this, m_nop_offsets[2], m_nops[2], m_nop_offsets[3], m_nops[3]){}

RdmIndsField &RdmIndsField::operator=(const RdmIndsField &other) {
    *this = static_cast<const base_t &>(other);
    return *this;
}

RdmIndsField &RdmIndsField::operator=(const FrmOnvConnection &conn) {
    m_frm = {conn.m_cre.inds(), conn.m_ann.inds()};
    return *this;
}

RdmIndsField &RdmIndsField::operator=(const BosOnvConnection &conn) {
    uint_t iind = 0ul;
    for (auto& pair : conn.m_cre.pairs()) {
        for (uint_t iop=0ul; iop<pair.m_nop; ++iop)
            m_bos.m_cre[iind++] = pair.m_imode;
    }
    iind = 0ul;
    for (auto& pair : conn.m_ann.pairs()) {
        for (uint_t iop=0ul; iop<pair.m_nop; ++iop)
            m_bos.m_ann[iind++] = pair.m_imode;
    }
    return *this;
}

RdmIndsField &RdmIndsField::operator=(const FrmBosOnvConnection &conn) {
    *this = conn.m_frm;
    *this = conn.m_bos;
    return *this;
}

bool RdmIndsField::is_ordered(const uint_t& iop, bool strict) const {
    return base_t::is_ordered(m_nop_offsets[iop], m_nop_offsets[iop]+m_nops[iop], strict, true);
}

bool RdmIndsField::is_ordered() const {
    return is_ordered(0, true) && is_ordered(1, true) && is_ordered(2, false) && is_ordered(3, false);
}

void RdmIndsField::common_frm_inds(uintv_t &common) const {
    common.clear();
    uint_t icre = 0ul;
    uint_t iann = 0ul;
    DEBUG_ASSERT_TRUE(is_ordered(), "indices are not properly ordered");
    while (icre < m_nops[0] && iann < m_nops[1]) {
        if (m_frm.m_cre[icre] > m_frm.m_ann[iann]) ++icre;
        else if (m_frm.m_cre[icre] < m_frm.m_ann[iann]) ++iann;
        else {
            // common element found
            common.push_back(m_frm.m_cre[icre]);
            ++icre; ++iann;
        }
    }
}
