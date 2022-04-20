//
// Created by rja on 07/04/2022.
//

#include "BasisData.h"

#include <utility>

FrmBasisData::FrmBasisData(size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved) :
        m_nsite(nsite), m_nspinorb(nsite*2), m_spin_resolved(spin_resolved), m_abgrp_map(std::move(abgrp_map)){}

FrmBasisData::FrmBasisData(size_t nsite, AbelianGroupMap abgrp_map) : FrmBasisData(nsite, abgrp_map, false){}

FrmBasisData::FrmBasisData(size_t nsite) : FrmBasisData(nsite, {nsite}, false){}

FrmBasisData::FrmBasisData() : FrmBasisData(0ul){}

bool FrmBasisData::operator==(const FrmBasisData &other) const {
    return m_nsite==other.m_nsite && m_abgrp_map==other.m_abgrp_map;
}


BosBasisData::BosBasisData(size_t nmode) : m_nmode(nmode){}

bool BosBasisData::operator==(const BosBasisData &other) const {
    return m_nmode==other.m_nmode;
}

BasisData::BasisData(FrmBasisData frm, BosBasisData bos) : m_frm(std::move(frm)), m_bos(std::move(bos)){}

void BasisData::require_pure_frm() const {
    REQUIRE_FALSE(m_bos.m_nmode, "MBF specification is not purely fermionic");
}

void BasisData::require_pure_bos() const {
    REQUIRE_FALSE(m_frm.m_nsite, "MBF specification is not purely bosonic");
}

bool BasisData::operator==(const BasisData &other) const {
    return m_frm==other.m_frm && m_bos==other.m_bos;
}
