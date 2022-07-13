//
// Created by Robert J. Anderson on 07/04/2022.
//

#include "BasisData.h"

#include <utility>


#if 0
//sys::frm::Size::Size(uint_t nsite)

uint_t sys::frm::Size::ncoeff_ind(bool restricted_orbs) const {
    return restricted_orbs ? m_nspinorb : m_nsite;
}

sys::Size::Size(uint_t nsite, uint_t nmode) : m_sites(nsite), m_nmode(nmode){}

void sys::Size::require_pure_frm() const {
    REQUIRE_FALSE(m_nmode, "Single particle basis specification is not purely fermionic");
}

void sys::Size::require_pure_bos() const {
    REQUIRE_FALSE(m_sites, "Single particle basis specification is not purely bosonic");
}

sys::frm::Basis::Basis(uint_t nelec, uint_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2) :
        m_nelec(nelec), m_sites(nsite), m_nvac(m_sites.m_nspinorb - nelec), m_spin_resolved(spin_resolved),
        m_abgrp_map(std::move(abgrp_map)), m_ms2(ms2),
        m_nelec_alpha(ms2_conserved() ? (m_nelec+m_ms2)/2 : 0ul),
        m_nelec_beta(ms2_conserved() ? m_nelec - m_nelec_alpha : 0ul),
        m_nvac_alpha(ms2_conserved() ? m_sites - m_nelec_alpha : 0ul),
        m_nvac_beta(ms2_conserved() ? m_sites - m_nelec_beta : 0ul){
    if (ms2_conserved())
        REQUIRE_EQ(uint_t(std::abs(m_ms2) % 2), m_nelec % 2, "2*Ms quantum number given incompatible with nelec");
}

sys::frm::Basis::Basis(uint_t nelec, uint_t nsite, int ms2) :
        sys::frm::Basis(nelec, nsite, {nsite}, false, ms2){}

sys::frm::Basis::Basis(uint_t nelec, uint_t nsite) :
        sys::frm::Basis(nelec, nsite, ~0){}

sys::frm::Basis::Basis(uint_t nsite) :
        sys::frm::Basis(0ul, nsite){}

sys::frm::Basis::Basis() : sys::frm::Basis(0ul){}

sys::frm::Basis::Basis(const sys::frm::Basis &hs1, const sys::frm::Basis &hs2) :
        sys::frm::Basis(hs1.m_nelec ? hs1.m_nelec : hs2.m_nelec,
                        hs1.m_sites ? hs1.m_sites : hs2.m_sites,
                        hs1.m_abgrp_map ? hs1.m_abgrp_map : hs2.m_abgrp_map,
                        hs1.m_spin_resolved || hs2.m_spin_resolved,
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

bool sys::frm::Basis::operator==(const sys::frm::Basis &other) const {
    return m_sites==other.m_sites && m_abgrp_map==other.m_abgrp_map;
}

BosHilbertSpace::BosHilbertSpace(uint_t nboson, uint_t nmode, bool nboson_conserve, uint_t occ_cutoff) :
        m_nboson(nboson), m_nmode(nmode), m_nboson_conserve(nboson_conserve), m_occ_cutoff(occ_cutoff){}

BosHilbertSpace::BosHilbertSpace(uint_t nboson, uint_t nmode) : BosHilbertSpace(nboson, nmode, true, max_bos_occ){}

BosHilbertSpace::BosHilbertSpace(uint_t nmode) : BosHilbertSpace(0ul, nmode, false, max_bos_occ){}

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

HilbertSpace::HilbertSpace(sys::frm::Basis frm, BosHilbertSpace bos) : m_frm(frm), m_bos(bos){}

HilbertSpace::HilbertSpace() : HilbertSpace(sys::frm::Basis(), BosHilbertSpace()){}

HilbertSpace::HilbertSpace(const HilbertSpace &hs1, const HilbertSpace &hs2) :
        HilbertSpace(sys::frm::Basis(hs1.m_frm, hs2.m_frm), BosHilbertSpace(hs1.m_bos, hs2.m_bos)){}


bool HilbertSpace::operator==(const HilbertSpace &other) const {
    return m_frm==other.m_frm && m_bos==other.m_bos;
}

sys::Size HilbertSpace::extents() const {
    return {m_frm.m_sites, m_bos.m_nmode};
}

void HilbertSpace::require_pure_frm() const {
    extents().require_pure_frm();
}

void HilbertSpace::require_pure_bos() const {
    extents().require_pure_bos();
}

HilbertSpace::operator bool() const {
    return m_frm || m_bos;
}

#endif

sys::frm::Size::Size(uint_t nsite) :
        m_nsite(nsite), m_nspinorb(2 * nsite), m_nspinorb_pair(integer::nspair(m_nspinorb)){}

uint_t sys::frm::Size::ncoeff_ind(bool spin_resolved) const {
    return spin_resolved ? m_nspinorb : m_nsite;
}

uint_t sys::frm::Size::ncoeff_ind(bool spin_resolved, uint_t nsite) {
    return spin_resolved ? 2*nsite : nsite;
}

sys::frm::Basis::Basis(uint_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved,
                       std::shared_ptr<lattice::Lattice> lattice) :
        Size(nsite), m_abgrp_map(std::move(abgrp_map)), m_spin_resolved(spin_resolved), m_lattice(lattice){
    if (*m_lattice) {
        REQUIRE_EQ(m_nsite, m_lattice->m_nsite, "incompatible lattice and site number");
    }
}

sys::frm::Basis::Basis(std::shared_ptr<lattice::Lattice> lattice) :
        Basis(lattice->m_nsite, {lattice->m_nsite}, false, lattice){
    REQUIRE_TRUE(*m_lattice, "cannot build a basis from a null lattice");
}

sys::frm::Basis::Basis(uint_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved) :
        Basis(nsite, abgrp_map, spin_resolved, lattice::make()){}

bool sys::frm::Basis::operator==(const sys::frm::Basis &other) const {
    return (m_nsite == other.m_nsite) && (m_abgrp_map == other.m_abgrp_map) && (m_spin_resolved == other.m_spin_resolved);
}

uint_t sys::frm::Basis::ncoeff_ind() const {
    return Size::ncoeff_ind(m_spin_resolved);
}

strmap_t sys::frm::Basis::info() const {
    strmap_t map;
    map.insert({"nsite", convert::to_string(m_nsite)});
    map.insert({"point group irreps", convert::to_string(m_abgrp_map.m_site_irreps)});
    map.insert({"spin resolved", convert::to_string(m_spin_resolved)});
    map.insert({"lattice", convert::to_string(m_lattice->m_info)});
    return map;
}

str_t sys::frm::Basis::to_string() const {
    return convert::to_string(info());
}

int sys::frm::Ms2::lowest_value(uint_t nelec) {
    return nelec&1ul;
}

sys::frm::Ms2::Ms2(int v, bool conserve) : conservation::Optional<int>(v, conserve, "2*Ms"){}

sys::frm::Ms2::Ms2() : conservation::Optional<int>("2*Ms"){}

sys::frm::Electrons::Electrons(uint_t n, sys::frm::Ms2 ms2) : m_n(n), m_npair(integer::nspair(m_n)), m_ms2(ms2),
                                                              m_nalpha(m_ms2.conserve() ? (m_n+m_ms2)/2 : 0ul),
                                                              m_nbeta(m_ms2.conserve() ? m_n-m_nalpha : 0ul) {
    if (m_ms2.conserve() && m_n)
        REQUIRE_EQ(uint_t(std::abs(m_ms2) % 2), m_n % 2,
                   "2*Ms quantum number given incompatible with number of electrons");
}

sys::frm::Electrons::Electrons(const sys::frm::Electrons &e1, const sys::frm::Electrons &e2) :
        Electrons(e1.m_n ? e1.m_n : e2.m_n, Ms2(e1.m_ms2, e2.m_ms2)){
    if (e1 && e2) REQUIRE_EQ(e1, e2, "incompatible numbers of electrons");
}

bool sys::frm::Electrons::operator==(const sys::frm::Electrons &other) const {
    return m_n==other.m_n && m_ms2==other.m_ms2;
}

std::map<str_t, str_t> sys::frm::Electrons::info() const {
    std::map<str_t, str_t> map;
    map.insert({"number", convert::to_string(m_n)});
    map.insert({"ms2", convert::to_string(m_ms2)});
    return map;
}

str_t sys::frm::Electrons::to_string() const {
    return convert::to_string(info());
}

sys::frm::Sector::Sector(sys::frm::Basis basis, sys::frm::Electrons elecs) :
        m_basis(basis), m_elecs(elecs),
        m_nvac(m_basis.m_nspinorb - m_elecs),
        m_nvac_alpha(m_elecs.m_ms2.conserve() ? m_basis.m_nsite - m_elecs.m_nalpha : 0ul),
        m_nvac_beta(m_elecs.m_ms2.conserve() ? m_basis.m_nsite - m_elecs.m_nbeta : 0ul){
    REQUIRE_LE_ALL(m_elecs, m_basis.m_nspinorb, "unphysical number of electrons");
}

bool sys::frm::Sector::operator==(const sys::frm::Sector &other) const {
    return m_basis==other.m_basis && m_elecs==other.m_elecs;
}

uint_t sys::frm::Sector::size() const {
    if (m_elecs.m_ms2.conserve()) {
        const auto na = integer::combinatorial(m_basis.m_nsite, m_elecs.m_nalpha);
        const auto nb = integer::combinatorial(m_basis.m_nsite, m_elecs.m_nbeta);
        return na * nb;
    }
    return integer::combinatorial(m_basis.m_nspinorb, m_elecs);
}

sys::bos::Size::Size(uint_t nmode) : m_nmode(nmode){}

sys::bos::Basis::Basis(uint_t nmode, uint_t occ_cutoff) : Size(nmode), m_occ_cutoff(occ_cutoff){}

bool sys::bos::Basis::operator==(const sys::bos::Basis &other) const {
    return (m_occ_cutoff == other.m_occ_cutoff) && (m_nmode == other.m_nmode);
}

strmap_t sys::bos::Basis::info() const {
    return {
            {"nmode",      convert::to_string(m_nmode)},
            {"occ_cutoff", convert::to_string(m_occ_cutoff)}
    };
}

str_t sys::bos::Basis::to_string() const {
    return convert::to_string(info());
}

sys::bos::Bosons::Bosons(uint_t v, bool conserve) : conservation::Optional<uint_t>(v, conserve, "nboson"){}

sys::bos::Bosons::Bosons() : Bosons(~0ul, true){}

sys::bos::Bosons::Bosons(const sys::bos::Bosons &b1, const sys::bos::Bosons &b2) : Bosons(b1.defined() ? b1 : b2){}

sys::bos::Sector::Sector(sys::bos::Basis basis, sys::bos::Bosons bosons) : m_basis(basis), m_bosons(bosons){}

bool sys::bos::Sector::operator==(const sys::bos::Sector &other) const {
    return m_basis==other.m_basis && m_bosons==other.m_bosons;
}

sys::Size::Size(uint_t nsite, uint_t nmode) : m_frm(nsite), m_bos(nmode){}

void sys::Size::require_pure_frm() const {
    REQUIRE_FALSE(m_bos, "Single particle basis specification is not purely fermionic");
}

void sys::Size::require_pure_bos() const {
    REQUIRE_FALSE(m_frm, "Single particle basis specification is not purely bosonic");
}

sys::Basis::Basis(sys::frm::Basis frm, sys::bos::Basis bos) : m_frm(std::move(frm)), m_bos(bos){}

void sys::Basis::require_pure_frm() const {
    size().require_pure_frm();
}

void sys::Basis::require_pure_bos() const {
    size().require_pure_bos();
}

sys::Size sys::Basis::size() const {
    return {static_cast<const frm::Size&>(m_frm), static_cast<const bos::Size&>(m_bos)};
}

sys::Sector::Sector(sys::frm::Sector frm, sys::bos::Sector bos) : m_frm(std::move(frm)), m_bos(std::move(bos)){}

sys::Sector::Sector(sys::frm::Sector frm) : Sector(frm, bos::Sector({0ul}, {})){}

sys::Sector::Sector(sys::bos::Sector bos) : Sector(frm::Sector({0ul}, {0ul}), bos){}

sys::Sector::Sector(sys::Basis basis, sys::Particles particles) :
        Sector(frm::Sector(basis.m_frm, particles.m_frm), bos::Sector(basis.m_bos, particles.m_bos)){}

sys::Basis sys::Sector::basis() const {
    return {m_frm.m_basis, m_bos.m_basis};
}

sys::Size sys::Sector::size() const {
    return basis().size();
}

sys::Particles sys::Sector::particles() const {
    return {m_frm.m_elecs, m_bos.m_bosons};
}
