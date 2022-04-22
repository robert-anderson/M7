//
// Created by rja on 07/04/2022.
//

#include "BasisData.h"

#include <utility>

FrmSites::FrmSites(size_t nsite) : m_nsite(nsite), m_nspinorb(2 * nsite){}

size_t FrmSites::ncoeff_ind(bool spin_resolved) const {
    return spin_resolved ? m_nspinorb : m_nsite;
}

BasisExtents::BasisExtents(size_t nsite, size_t nmode) : m_sites(nsite), m_nmode(nmode){}

void BasisExtents::require_pure_frm() const {
    REQUIRE_FALSE(m_nmode, "Single particle basis specification is not purely fermionic");
}

void BasisExtents::require_pure_bos() const {
    REQUIRE_FALSE(m_sites, "Single particle basis specification is not purely bosonic");
}

FrmHilbertSpace::FrmHilbertSpace(size_t nelec, size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2,
                                 bool ms2_conserve) :
        m_nelec(nelec), m_sites(nsite), m_spin_resolved(spin_resolved),
        m_abgrp_map(std::move(abgrp_map)), m_ms2(ms2), m_ms2_conserve(ms2_conserve){}

FrmHilbertSpace::FrmHilbertSpace(size_t nelec, size_t nsite) :
        FrmHilbertSpace(nelec, nsite, {nsite}, false, 0, false){}

bool FrmHilbertSpace::operator==(const FrmHilbertSpace &other) const {
    return m_sites==other.m_sites && m_abgrp_map==other.m_abgrp_map;
}

BosHilbertSpace::BosHilbertSpace(size_t nmode, size_t nboson, bool nboson_conserve, size_t occ_cutoff) :
        m_nmode(nmode), m_nboson(nboson), m_nboson_conserve(nboson_conserve), m_occ_cutoff(occ_cutoff){}

BosHilbertSpace::BosHilbertSpace(size_t nmode) : BosHilbertSpace(nmode, 0ul, false, defs::max_bos_occ){}

bool BosHilbertSpace::operator==(const BosHilbertSpace &other) const {
    return m_nmode==other.m_nmode;
}

HilbertSpace::HilbertSpace(FrmHilbertSpace frm, BosHilbertSpace bos) :
        m_frm(frm), m_bos(bos), m_extents(frm.m_sites, bos.m_nmode){}

bool HilbertSpace::operator==(const HilbertSpace &other) const {
    return m_frm==other.m_frm && m_bos==other.m_bos;
}

void HilbertSpace::require_pure_frm() const {
    m_extents.require_pure_frm();
}

void HilbertSpace::require_pure_bos() const {
    m_extents.require_pure_bos();
}
