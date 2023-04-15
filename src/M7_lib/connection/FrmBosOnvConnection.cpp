//
// Created by Robert J. Anderson on 10/08/2021.
//

#include "FrmBosOnvConnection.h"
#include "ComOps.h"

FrmBosOnvConnection::FrmBosOnvConnection(sys::Size size) :
    m_frm(size.m_frm), m_bos(size.m_bos){}

FrmBosOnvConnection::FrmBosOnvConnection(const FrmBosOnvField& mbf) :
        FrmBosOnvConnection({mbf.m_frm.m_basis.m_nsite, mbf.m_bos.m_basis.m_nmode}){}


void FrmBosOnvConnection::clear() {
    m_frm.clear();
    m_bos.clear();
}

void FrmBosOnvConnection::connect(const FrmOnvField& src, const FrmOnvField& dst) {
    m_bos.clear();
    m_frm.connect(src, dst);
}

bool FrmBosOnvConnection::connect(const FrmOnvField& src, const FrmOnvField& dst, com_ops::Frm& com) {
    m_bos.clear();
    return m_frm.connect(src, dst, com);
}

void FrmBosOnvConnection::connect(const BosOnvField& src, const BosOnvField& dst) {
    m_frm.clear();
    m_bos.connect(src, dst);
}

bool FrmBosOnvConnection::connect(const BosOnvField& src, const BosOnvField& dst, com_ops::Bos& com) {
    m_frm.clear();
    return m_bos.connect(src, dst, com);
}

void FrmBosOnvConnection::connect(const FrmBosOnvField& src, const FrmBosOnvField& dst) {
    m_frm.connect(src.m_frm, dst.m_frm);
    m_bos.connect(src.m_bos, dst.m_bos);
}

bool FrmBosOnvConnection::connect(const FrmBosOnvField& src, const FrmBosOnvField& dst, com_ops::FrmBos& com) {
    m_bos.connect(src.m_bos, dst.m_bos);
    return m_frm.connect(src.m_frm, dst.m_frm, com.m_frm);
}


void FrmBosOnvConnection::apply(const FrmBosOnvField& src, FrmBosOnvField& dst) const {
    m_frm.apply(src.m_frm, dst.m_frm);
    m_bos.apply(src.m_bos, dst.m_bos);
}

OpSig FrmBosOnvConnection::exsig() const {
    return {{m_frm.m_cre.size(), m_frm.m_ann.size()}, {m_bos.m_cre.size(), m_bos.m_ann.size()}};
}

bool FrmBosOnvConnection::respects_occ_range(const FrmBosOnvField& src, uint_t nboson_max) const {
    return m_bos.respects_occ_range(src.m_bos, nboson_max);
}
