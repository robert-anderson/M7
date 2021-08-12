//
// Created by rja on 11/08/2021.
//

#include "MaeTable.h"

std::array<size_t, 4> MaeInds::make_nops() const {
    std::array<size_t, 4> nops{};
    nops[0] = conn_utils::extract_ncref(m_exsig);
    nops[1] = conn_utils::extract_nannf(m_exsig);
    nops[2] = conn_utils::extract_ncreb(m_exsig);
    nops[3] = conn_utils::extract_nannb(m_exsig);
    return nops;
}

std::array<size_t, 4> MaeInds::make_nop_offsets() const {
    auto nop_offsets = m_nops;
    for (size_t i = 1ul; i < 4; ++i) nop_offsets[i] = nop_offsets[i - 1] + m_nops[i - 1];
    return nop_offsets;
}

MaeInds::MaeInds(Row *row, size_t exsig, std::string name) :
        fields::Numbers<defs::mev_ind_t, 1>(row, {conn_utils::extract_nop(exsig)}, name),
        m_exsig(exsig), m_nops(make_nops()), m_nop_offsets(make_nop_offsets()) {}

MaeInds &MaeInds::operator=(const std::array<defs::inds, 4> &inds) {
    for (size_t iop = 0ul; iop<4; ++iop) {
        const size_t offset = m_nop_offsets[iop];
        for (size_t iind=0ul; iind<inds[iop].size(); ++iind) (*this)[offset+iind] = inds[iop][iind];
    }
    return *this;
}

const defs::mev_ind_t &MaeInds::get_ref(const size_t &iop, const size_t &iind) const {
    DEBUG_ASSERT_LT(iop, 4ul, "operator type index OOB");
    DEBUG_ASSERT_LT(iind, m_nops[iop], "operator index OOB");
    return (*this)[m_nop_offsets[iop]+iind];
}

defs::mev_ind_t &MaeInds::get_ref(const size_t &iop, const size_t &iind) {
    DEBUG_ASSERT_LT(iop, 4ul, "operator type index OOB");
    DEBUG_ASSERT_LT(iind, m_nops[iop], "operator index OOB");
    return (*this)[m_nop_offsets[iop]+iind];
}

bool MaeInds::is_ordered(const size_t& iop, bool strict) const {
    return base_t::is_ordered(m_nop_offsets[iop], m_nop_offsets[iop]+m_nops[iop], strict, true);
}

bool MaeInds::is_ordered() const {
    return is_ordered(0, true) && is_ordered(1, true) && is_ordered(2, false) && is_ordered(3, false);
}

void MaeInds::common_frm_inds(defs::inds &common) const {
    common.clear();
    size_t icre = 0ul;
    size_t iann = 0ul;
    DEBUG_ASSERT_TRUE(is_ordered(), "indices are not properly ordered");
    while (icre < m_nops[0] && iann < m_nops[1]) {
        if (get_ref<0>(icre) > get_ref<1>(iann)) ++icre;
        else if (get_ref<0>(icre) < get_ref<1>(iann)) ++iann;
        else {
            // common element found
            common.push_back((*this)[icre]);
            ++icre; ++iann;
        }
    }
}

std::string MaeInds::get_exsig_string() const {
    return log::format("{}{}{}{}",
                       conn_utils::extract_ncref(m_exsig), conn_utils::extract_nannf(m_exsig),
                       conn_utils::extract_ncreb(m_exsig), conn_utils::extract_nannb(m_exsig));
}

const size_t &MaeInds::nop(const size_t &iop) const {
    DEBUG_ASSERT_LT(iop, 4ul, "operator type index OOB");
    return m_nops[iop];
}

SpecMomInds::SpecMomInds(Row *row, size_t exsig, std::string name) :
        base_t(row, {nullptr, exsig, "left " + name}, {nullptr, exsig, "right " + name}),
        m_exsig(exsig), m_name(name), m_left(get<0>()), m_right(get<1>()) {}

SpecMomInds::SpecMomInds(const SpecMomInds &other) : SpecMomInds(other.row_of_copy(), other.m_exsig, other.m_name) {}
