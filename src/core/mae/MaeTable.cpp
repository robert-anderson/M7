//
// Created by rja on 11/08/2021.
//

#include "MaeTable.h"

mae_inds::Frm::Frm(Row *row, size_t nann, size_t ncre, std::string name) :
        base_t(row, {nullptr, {nann}, name.empty() ? "" : name + "_ann"},
               {nullptr, {ncre}, name.empty() ? "" : name + "_cre"}),
        m_name(name), m_ann(get<0>()), m_cre(get<1>()) {
}

mae_inds::Frm::Frm(const Frm &other) :
        Frm(other.m_ann.row_of_copy(),
                       other.m_ann.m_format.m_nelement,
                       other.m_cre.m_format.m_nelement, other.m_name) {}

mae_inds::Frm &mae_inds::Frm::operator=(const Frm &other) {
    m_ann = other.m_ann;
    m_cre = other.m_cre;
    return *this;
}

mae_inds::Frm &mae_inds::Frm::operator=(const std::pair<defs::inds, defs::inds> &inds) {
    for (size_t i = 0ul; i < inds.first.size(); ++i) m_ann[i] = inds.first[i];
    for (size_t i = 0ul; i < inds.second.size(); ++i) m_cre[i] = inds.second[i];
    return *this;
}

bool mae_inds::Frm::is_ordered() const {
    return m_ann.is_ordered(false, true) && m_cre.is_ordered(false, true);
}

size_t mae_inds::Frm::ncommon_sq_op_ind() const {
    size_t ncommon = 0ul;
    size_t icre = 0ul;
    size_t iann = 0ul;
    ASSERT(is_ordered());
    while (icre < m_cre.nelement() && iann < m_ann.nelement()) {
        if (m_cre[icre] > m_ann[iann]) ++icre;
        else if (m_cre[icre] < m_ann[iann]) ++iann;
        else {
            // common element found
            ++ncommon;
            ++icre;
            ++iann;
        }
    }
    return ncommon;
}

void mae_inds::Frm::common_sq_op_inds(defs::inds &common) const {
    common.clear();
    size_t icre = 0ul;
    size_t iann = 0ul;
    ASSERT(is_ordered());
    while (icre < m_cre.nelement() && iann < m_ann.nelement()) {
        if (m_cre[icre] > m_ann[iann]) ++icre;
        else if (m_cre[icre] < m_ann[iann]) ++iann;
        else {
            // common element found
            common.push_back(m_cre[icre]);
            ++icre;
            ++iann;
        }
    }
}
