//
// Created by Robert J. Anderson on 24/04/22.
//

#include "Suites.h"
#include "BasisData.h"

suite::MbfsRow::MbfsRow(const sys::Basis &basis) :
        m_frm(this, basis.m_frm, "fermion ONV"),
        m_frmbos(this, basis, "fermion-boson ONV"),
        m_bos(this, basis.m_bos, "boson ONV"){}

suite::MbfsRow::MbfsRow(const sys::Sector &sector) : MbfsRow(sector.basis()){}


suite::Mbfs::Mbfs(const sys::Basis &basis) :
    buffered::Table<MbfsRow>(basis) {
    m_row.push_back_jump();
}

suite::Mbfs::Mbfs(const sys::Sector &sector) : Mbfs(sector.basis()){}

suite::Conns::Conns(const sys::Size &size) :
        m_frmonv(size.m_frm), m_bosonv(size.m_bos), m_frmbosonv(size){}

suite::ComOps::ComOps(const sys::Size &size) : m_frm(size.m_frm), m_frmbos(size), m_bos(size.m_bos){}