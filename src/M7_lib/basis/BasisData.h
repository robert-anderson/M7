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
        DEBUG_ASSERT_LT(ispin, 2, "spin channel index OOB");
        DEBUG_ASSERT_LT(isite, 2, "site index OOB");
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
     * @param spin_resolved
     *  true if the basis is spin-resolved e.g. UHF
     * @return
     *  number of indices needed in the access of coefficients
     */
    size_t ncoeff_ind(bool spin_resolved) const;
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
     * true if the two spin orbitals corresponding to the same site have identical functional form e.g. in UHF basis
     */
    const bool m_spin_resolved;
    /**
     * mapping from fermion site indices to Abelian group labels
     */
    const AbelianGroupMap m_abgrp_map;
    /**
     * the 2*Ms value to which the fermionic many-body basis is constrained
     */
    const int m_ms2;
    /**
     * if true, the above ms2 value is respected, else it is not enforced upon the MBFs, but still influences the choice
     * of initial states and references.
     */
    const bool m_ms2_conserve;

    FrmHilbertSpace(size_t nelec, size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2, bool ms2_conserve);
    /*
     * non-resolved spin, and no Ms2 restriction
     */
    FrmHilbertSpace(size_t nelec, size_t nsite);

    bool operator==(const FrmHilbertSpace& other) const;

    operator bool() const {
        return m_sites;
    }
};

/**
 * properties of the bosonic Hilbert space
 */
struct BosHilbertSpace {
    /**
     * number of bosonic modes
     */
    const size_t m_nmode;
    /**
     * number of bosons in a number-conserving system
     */
    const size_t m_nboson;
    /**
     * if true, the above nboson value is respected, else it is not enforced upon the MBFs, but still influences the
     * choice of initial states and references.
     */
    const bool m_nboson_conserve;
    /**
     * number of bosons permitted to occupy any given mode
     */
    const size_t m_occ_cutoff;

    BosHilbertSpace(size_t nmode, size_t nboson, bool nboson_conserve, size_t occ_cutoff);

    BosHilbertSpace(size_t nmode);

    bool operator==(const BosHilbertSpace& other) const;

    operator bool() const {
        return m_nmode;
    }
};

/**
 * admits a common interface for all initialisations of MBFs
 */
struct HilbertSpace {
    const FrmHilbertSpace m_frm;
    const BosHilbertSpace m_bos;
    const BasisExtents m_extents;

    HilbertSpace(FrmHilbertSpace frm, BosHilbertSpace bos);

    bool operator==(const HilbertSpace& other) const;

    void require_pure_frm() const;
    void require_pure_bos() const;

    operator bool() const {
        return m_frm || m_bos;
    }
};

#endif //M7_BASISDATA_H