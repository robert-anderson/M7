//
// Created by rja on 24/08/2021.
//

#include "CachedOrbs.h"

CachedOrbs::CachedOrbs(const AbelianGroupMap &grp_map) : m_occ(grp_map), m_vac(grp_map){
    clear();
}

void CachedOrbs::clear() {
    m_occ.clear();
    m_vac.clear();
    m_nonempty_pair_labels.clear();
}

const SpinSymOccOrbs &CachedOrbs::occ(const field::FrmOnv &mbf) {
    if (m_occ.empty()) m_occ.update(mbf);
    return m_occ;
}

const SpinSymVacOrbs &CachedOrbs::vac(const field::FrmOnv &mbf) {
    if (m_vac.empty()) m_vac.update(mbf);
    return m_vac;
}

const defs::inds &CachedOrbs::nonempty_pair_labels(const field::FrmOnv &mbf) {
    if (m_nonempty_pair_labels.empty()){
        auto& occ = this->occ(mbf);
        auto& vac = this->vac(mbf);
        for (size_t i=0ul; i<occ.m_format.m_nelement; ++i){
            if (!occ[i].empty() && !vac[i].empty()) m_nonempty_pair_labels.push_back(i);
        }
    }
    return m_nonempty_pair_labels;
}
