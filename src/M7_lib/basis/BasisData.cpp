//
// Created by rja on 07/04/2022.
//

#include "BasisData.h"

#include <utility>

FrmBasisData::FrmBasisData(size_t nsite, AbelianGroupMap abgrp_map) :
        m_nsite(nsite), m_nspinorb(nsite*2), m_abgrp_map(std::move(abgrp_map)){}

FrmBasisData::FrmBasisData(size_t nsite) : FrmBasisData(nsite, {nsite}){}

FrmBasisData::FrmBasisData() : FrmBasisData(0ul){}

bool FrmBasisData::operator==(const FrmBasisData &other) const {
    return m_nsite==other.m_nsite && m_abgrp_map==other.m_abgrp_map;
}

BosBasisData::BosBasisData(size_t nmode, size_t nboson_max) : m_nmode(nmode), m_nboson_max(nboson_max){
    REQUIRE_LE(nboson_max, defs::max_bos_occ,
               "the given nboson_max exceeds the limits of the defs::bos_occ_t container");
}

bool BosBasisData::operator==(const BosBasisData &other) const {
    return m_nmode==other.m_nmode && m_nboson_max==other.m_nboson_max;
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
