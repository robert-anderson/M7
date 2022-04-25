//
// Created by rja on 10/08/2021.
//

#include "FrmBosOnvConnection.h"
#include "ComOps.h"

FrmBosOnvConnection::FrmBosOnvConnection(sys::Size extents) :
    m_frm(extents.m_sites), m_bos(extents.m_nmode){}

FrmBosOnvConnection::FrmBosOnvConnection(const FrmBosOnvField &mbf) :
        FrmBosOnvConnection({mbf.m_frm.m_hs.m_sites, mbf.m_bos.m_hs.m_nmode}){}


void FrmBosOnvConnection::clear() {
    m_frm.clear();
    m_bos.clear();
}

void FrmBosOnvConnection::connect(const FrmBosOnvField &src, const FrmBosOnvField &dst) {
    m_frm.connect(src.m_frm, dst.m_frm);
    m_bos.connect(src.m_bos, dst.m_bos);
}

bool FrmBosOnvConnection::connect(const FrmBosOnvField &src, const FrmBosOnvField &dst, com_ops::FrmBos &com) {
    m_bos.connect(src.m_bos, dst.m_bos);
    return m_frm.connect(src.m_frm, dst.m_frm, com.m_frm);
}


void FrmBosOnvConnection::apply(const FrmBosOnvField &src, FrmBosOnvField &dst) const {
    m_frm.apply(src.m_frm, dst.m_frm);
    m_bos.apply(src.m_bos, dst.m_bos);
}

size_t FrmBosOnvConnection::exsig() const {
    return exsig_utils::encode(m_frm.m_cre.size(), m_frm.m_ann.size(), m_bos.m_cre.size(), m_bos.m_ann.size());
}

bool FrmBosOnvConnection::respects_occ_range(const FrmBosOnvField &src, size_t nboson_max) const {
    return m_bos.respects_occ_range(src.m_bos, nboson_max);
}
