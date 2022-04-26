//
// Created by rja on 24/04/22.
//

#include "Suites.h"
#include "BasisData.h"

suite::MbfsRow::MbfsRow(const sys::Sector &sector) :
        m_frm(this, sector, "fermion ONV"),
        m_frmbos(this, sector, "fermion-boson ONV"),
        m_bos(this, sector, "boson ONV"){}

suite::Mbfs::Mbfs(const sys::Sector &sector) : BufferedTable<MbfsRow>("Work space for MBFs", {{sector}}){
    m_row.push_back_jump();
}

suite::Conns::Conns(const sys::Size &size) :
        m_frmonv(size.m_frm), m_bosonv(size.m_bos), m_frmbosonv(size){}

suite::ComOps::ComOps(const sys::Size &size) : m_frm(size.m_frm), m_frmbos(size), m_bos(size.m_bos){}