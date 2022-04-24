//
// Created by rja on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"


class FrmSites {
    const size_t m_nsite;
public:
    const size_t m_nspinorb;
    FrmSites(size_t nsite);

    operator const size_t& () const {
        return m_nsite;
    }

    size_t isite(size_t ispinorb) const {
        DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
        return ispinorb < m_nsite ? ispinorb : ispinorb-m_nsite;
    }
    size_t ispin(size_t ispinorb) const {
        DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
        return ispinorb >= m_nsite;
    }
    size_t ispinorb(size_t ispin, size_t isite) const {
        DEBUG_ASSERT_LT(ispin, 2ul, "spin channel index OOB");
        DEBUG_ASSERT_LT(isite, m_nsite, "site index OOB");
        return ispin ? isite + m_nsite : isite;
    }
    size_t ispinorb(std::pair<size_t, size_t> pair) const {
        return ispinorb(pair.first, pair.second);
    }
    int ms2(size_t ispinorb) const {
        DEBUG_ASSERT_LT(ispinorb, m_nspinorb, "spin orbital index OOB");
        return ispinorb < m_nsite ? 1 : -1;
    }
    /**
     * @param restricted_orbs
     *  true if the basis is spin-resolved e.g. UHF
     * @return
     *  number of indices needed in the access of coefficients
     */
    size_t ncoeff_ind(bool restricted_orbs) const;
};

struct BasisExtents {
    const FrmSites m_sites;
    const size_t m_nmode;
    BasisExtents(size_t nsite, size_t nmode);
    void require_pure_frm() const;
    void require_pure_bos() const;
};

/**
 * properties of the fermionic Hilbert space
 */
struct FrmHilbertSpace {
    /**
     * number of electrons in the system
     */
    const size_t m_nelec;
    /**
     * number of sites and spin orbitals
     */
    const FrmSites m_sites;
    /**
     * number of vacant spin orbitals
     */
    const size_t m_nvac;
    /**
     * true if the two spin orbitals corresponding to the same site have identical functional form e.g. in UHF basis
     */
    const bool m_restricted_orbs;
    /**
     * mapping from fermion site indices to Abelian group labels
     */
    const AbelianGroupMap m_abgrp_map;
    /**
     * the 2*Ms value to which the fermionic many-body basis is constrained (~0 if conserved, i.e. commutator [H, Sz]!=0)
     */
    const int m_ms2;
    /**
     * numbers of occupied alpha and beta electrons (0 if spin unconserved)
     */
    const size_t m_nelec_alpha, m_nelec_beta;
    /**
     * numbers of occupied alpha and beta vacant orbitals (0 if spin unconserved)
     */
    const size_t m_nvac_alpha, m_nvac_beta;

    FrmHilbertSpace(size_t nelec, size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2);
    /*
     * non-resolved spin, C1 point group (no spatial symmetry)
     */
    FrmHilbertSpace(size_t nelec, size_t nsite, int ms2);
    /*
     * non-resolved spin, C1 point group (no spatial symmetry), and undefined Ms2 restriction
     */
    FrmHilbertSpace(size_t nelec, size_t nsite);
    FrmHilbertSpace();

    /*
     * combine the properties of two Hilbert spaces
     */
    FrmHilbertSpace(const FrmHilbertSpace& hs1, const FrmHilbertSpace& hs2);

    bool operator==(const FrmHilbertSpace& other) const;

    operator bool() const {
        return m_sites;
    }

    bool ms2_conserved() const {
        return m_ms2!=~0;
    }
};

/**
 * properties of the bosonic Hilbert space
 */
struct BosHilbertSpace {
    /**
     * number of bosons in a number-conserving system
     */
    const size_t m_nboson;
    /**
     * number of bosonic modes
     */
    const size_t m_nmode;
    /**
     * if true, the above nboson value is respected, else it is not enforced upon the MBFs, but still influences the
     * choice of initial states and references.
     */
    const bool m_nboson_conserve;
    /**
     * number of bosons permitted to occupy any given mode
     */
    const size_t m_occ_cutoff;

    BosHilbertSpace(size_t nboson, size_t nmode, bool nboson_conserve, size_t occ_cutoff);
    /*
     * number non-conserved, maximum occupation cutoff
     */
    BosHilbertSpace(size_t nmode);
    BosHilbertSpace();
    /*
     * combine the properties of two Hilbert spaces
     */
    BosHilbertSpace(const BosHilbertSpace& hs1, const BosHilbertSpace& hs2);

    bool operator==(const BosHilbertSpace& other) const;

    operator bool() const;
};

/**
 * admits a common interface for all initialisations of MBFs
 */
struct HilbertSpace {
    const FrmHilbertSpace m_frm;
    const BosHilbertSpace m_bos;
    const BasisExtents m_extents;

    HilbertSpace(FrmHilbertSpace frm, BosHilbertSpace bos);
    HilbertSpace();
    /*
     * combine the properties of two Hilbert spaces
     */
    HilbertSpace(const HilbertSpace& hs1, const HilbertSpace& hs2);

    bool operator==(const HilbertSpace& other) const;

    void require_pure_frm() const;
    void require_pure_bos() const;

    operator bool() const;
};

#endif //M7_BASISDATA_H