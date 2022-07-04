//
// Created by Robert J. Anderson on 12/08/2021.
//

#include "MaeIndsField.h"
#include <M7_lib/util/Exsig.h>


uinta_t<4> MaeIndsField::make_nops() const {
    uinta_t<4> nops{};
    nops[0] = exsig::decode_nfrm_cre(m_exsig);
    nops[1] = exsig::decode_nfrm_ann(m_exsig);
    nops[2] = exsig::decode_nbos_cre(m_exsig);
    nops[3] = exsig::decode_nbos_ann(m_exsig);
    return nops;
}

uinta_t<4> MaeIndsField::make_nop_offsets() const {
    uinta_t<4> nop_offsets{};
    for (uint_t i = 1ul; i < 4; ++i) nop_offsets[i] = nop_offsets[i - 1] + m_nops[i - 1];
    return nop_offsets;
}

MaeIndsField::MaeIndsField(Row *row, uint_t exsig, str_t name) :
    NdNumberField<mev_ind_t, 1>(row, {exsig::decode_nop(exsig)}, name),
    m_exsig(exsig), m_nops(make_nops()), m_nop_offsets(make_nop_offsets()),
    m_frm(*this, m_nop_offsets[0], m_nops[0], m_nop_offsets[1], m_nops[1]),
    m_bos(*this, m_nop_offsets[2], m_nops[2], m_nop_offsets[3], m_nops[3]){}

MaeIndsField::MaeIndsField(const MaeIndsField &other) :
    NdNumberField<mev_ind_t, 1>(other), m_exsig(other.m_exsig), m_nops(make_nops()), m_nop_offsets(make_nop_offsets()),
    m_frm(*this, m_nop_offsets[0], m_nops[0], m_nop_offsets[1], m_nops[1]),
    m_bos(*this, m_nop_offsets[2], m_nops[2], m_nop_offsets[3], m_nops[3]){}

MaeIndsField &MaeIndsField::operator=(const MaeIndsField &other) {
    *this = static_cast<const base_t &>(other);
    return *this;
}

MaeIndsField &MaeIndsField::operator=(const FrmOnvConnection &conn) {
    m_frm = {conn.m_cre.inds(), conn.m_ann.inds()};
    return *this;
}

MaeIndsField &MaeIndsField::operator=(const BosOnvConnection &conn) {
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

MaeIndsField &MaeIndsField::operator=(const FrmBosOnvConnection &conn) {
    *this = conn.m_frm;
    *this = conn.m_bos;
    return *this;
}

bool MaeIndsField::is_ordered(const uint_t& iop, bool strict) const {
    return base_t::is_ordered(m_nop_offsets[iop], m_nop_offsets[iop]+m_nops[iop], strict, true);
}

bool MaeIndsField::is_ordered() const {
    return is_ordered(0, true) && is_ordered(1, true) && is_ordered(2, false) && is_ordered(3, false);
}

void MaeIndsField::common_frm_inds(uintv_t &common) const {
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
