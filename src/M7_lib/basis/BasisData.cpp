//
// Created by rja on 07/04/2022.
//

#include "BasisData.h"

#include <utility>

FrmSites::FrmSites(size_t nsite) : m_nsite(nsite), m_nspinorb(2 * nsite){}

size_t FrmSites::ncoeff_ind(bool restricted_orbs) const {
    return restricted_orbs ? m_nspinorb : m_nsite;
}

BasisExtents::BasisExtents(size_t nsite, size_t nmode) : m_sites(nsite), m_nmode(nmode){}

void BasisExtents::require_pure_frm() const {
    REQUIRE_FALSE(m_nmode, "Single particle basis specification is not purely fermionic");
}

void BasisExtents::require_pure_bos() const {
    REQUIRE_FALSE(m_sites, "Single particle basis specification is not purely bosonic");
}

FrmHilbertSpace::FrmHilbertSpace(size_t nelec, size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2) :
        m_nelec(nelec), m_sites(nsite), m_nvac(m_sites.m_nspinorb - nelec), m_restricted_orbs(spin_resolved),
        m_abgrp_map(std::move(abgrp_map)), m_ms2(ms2),
        m_nelec_alpha(ms2_conserved() ? (m_nelec+m_ms2)/2 : 0ul),
        m_nelec_beta(ms2_conserved() ? m_nelec - m_nelec_alpha : 0ul),
        m_nvac_alpha(ms2_conserved() ? m_sites - m_nelec_alpha : 0ul),
        m_nvac_beta(ms2_conserved() ? m_sites - m_nelec_beta : 0ul){
    if (ms2_conserved())
        REQUIRE_EQ(size_t(std::abs(m_ms2) % 2), m_nelec % 2, "2*Ms quantum number given incompatible with nelec");
}

FrmHilbertSpace::FrmHilbertSpace(size_t nelec, size_t nsite, int ms2) :
        FrmHilbertSpace(nelec, nsite, {nsite}, false, ms2){}

FrmHilbertSpace::FrmHilbertSpace(size_t nelec, size_t nsite) :
        FrmHilbertSpace(nelec, nsite, ~0){}

FrmHilbertSpace::FrmHilbertSpace(size_t nsite) :
        FrmHilbertSpace(0ul, nsite){}

FrmHilbertSpace::FrmHilbertSpace() : FrmHilbertSpace(0ul){}

FrmHilbertSpace::FrmHilbertSpace(const FrmHilbertSpace &hs1, const FrmHilbertSpace &hs2) :
        FrmHilbertSpace(hs1.m_nelec ? hs1.m_nelec : hs2.m_nelec,
                        hs1.m_sites ? hs1.m_sites : hs2.m_sites,
                        hs1.m_abgrp_map ? hs1.m_abgrp_map : hs2.m_abgrp_map,
                        hs1.m_restricted_orbs || hs2.m_restricted_orbs,
                        hs1.m_ms2 != ~0 ? hs1.m_ms2 : hs2.m_ms2) {
    if (hs1 && hs2) {
        /*
         * both Hilbert spaces are non-null so compatibility needs to be checked
         */
        REQUIRE_EQ(hs1.m_nelec, hs2.m_nelec, "incompatible numbers of electrons");
        REQUIRE_EQ(hs1.m_sites, hs2.m_sites, "incompatible numbers of sites");
        if (hs1.m_abgrp_map && hs2.m_abgrp_map){
            // both Hilbert spaces define point groups, so check compatibility
            REQUIRE_EQ(hs1.m_abgrp_map, hs2.m_abgrp_map, "incompatible Abelian group maps");
        }
        if (hs1.m_ms2 != ~0 && hs2.m_ms2 != ~0){
            // both Hilbert spaces define spin sector restrictions, so check compatibility
            REQUIRE_EQ(hs1.m_ms2, hs2.m_ms2, "incompatible spin sector restrictions");
        }
    }
}

bool FrmHilbertSpace::operator==(const FrmHilbertSpace &other) const {
    return m_sites==other.m_sites && m_abgrp_map==other.m_abgrp_map;
}

BosHilbertSpace::BosHilbertSpace(size_t nboson, size_t nmode, bool nboson_conserve, size_t occ_cutoff) :
        m_nboson(nboson), m_nmode(nmode), m_nboson_conserve(nboson_conserve), m_occ_cutoff(occ_cutoff){}

BosHilbertSpace::BosHilbertSpace(size_t nboson, size_t nmode) : BosHilbertSpace(nboson, nmode, true, defs::max_bos_occ){}

BosHilbertSpace::BosHilbertSpace(size_t nmode) : BosHilbertSpace(0ul, nmode, false, defs::max_bos_occ){}

BosHilbertSpace::BosHilbertSpace() : BosHilbertSpace(0ul){}

BosHilbertSpace::BosHilbertSpace(const BosHilbertSpace &hs1, const BosHilbertSpace &hs2) :
        BosHilbertSpace(
                hs1.m_nboson ? hs1.m_nboson : hs2.m_nboson,
                hs1.m_nmode ? hs1.m_nmode : hs2.m_nmode,
                hs1.m_nboson_conserve && hs2.m_nboson_conserve,
                hs1.m_occ_cutoff ? hs1.m_occ_cutoff : hs2.m_occ_cutoff) {
    if (hs1 && hs2){
        /*
         * both Hilbert spaces are non-null so compatibility needs to be checked
         */
        REQUIRE_EQ(hs1.m_nboson, hs2.m_nboson, "incompatible numbers of bosons");
        REQUIRE_EQ(hs1.m_nmode, hs2.m_nmode, "incompatible numbers of boson modes");
        REQUIRE_EQ(hs1.m_occ_cutoff, hs2.m_occ_cutoff, "incompatible mode occupancy cutoffs");
    }
}

bool BosHilbertSpace::operator==(const BosHilbertSpace &other) const {
    return m_nmode==other.m_nmode;
}

BosHilbertSpace::operator bool() const {
    return m_nmode;
}

HilbertSpace::HilbertSpace(FrmHilbertSpace frm, BosHilbertSpace bos) :
        m_frm(frm), m_bos(bos), m_extents(frm.m_sites, bos.m_nmode){}

HilbertSpace::HilbertSpace() : HilbertSpace(FrmHilbertSpace(), BosHilbertSpace()){}

HilbertSpace::HilbertSpace(const HilbertSpace &hs1, const HilbertSpace &hs2) :
        HilbertSpace(FrmHilbertSpace(hs1.m_frm, hs2.m_frm), BosHilbertSpace(hs1.m_bos, hs2.m_bos)){}


bool HilbertSpace::operator==(const HilbertSpace &other) const {
    return m_frm==other.m_frm && m_bos==other.m_bos;
}

void HilbertSpace::require_pure_frm() const {
    m_extents.require_pure_frm();
}

void HilbertSpace::require_pure_bos() const {
    m_extents.require_pure_bos();
}

HilbertSpace::operator bool() const {
    return m_frm || m_bos;
}