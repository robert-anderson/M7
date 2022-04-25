//
// Created by rja on 24/04/22.
//

#include "Suites.h"

suite::MbfsRow::MbfsRow(const HilbertSpace &hs) :
        m_frm(this, hs, "fermion ONV"),
        m_frmbos(this, hs, "fermion-boson ONV"),
        m_bos(this, hs, "boson ONV"){}

suite::Mbfs::Mbfs(const HilbertSpace &hs) : BufferedTable<MbfsRow>("Work space for MBFs", {{hs}}){
    m_row.push_back_jump();
}

suite::Conns::Conns(const BasisExtents &extents) :
        m_frmonv(extents.m_sites), m_bosonv(extents.m_nmode), m_frmbosonv(extents){}

suite::Conns::Conns(const HilbertSpace &hs) : Conns(hs.extents()){}

suite::ComOps::ComOps(const BasisExtents &extents) : m_frm(extents.m_sites), m_frmbos(extents), m_bos(extents.m_nmode){}

suite::ComOps::ComOps(const HilbertSpace &hs) : ComOps(hs.extents()){}
