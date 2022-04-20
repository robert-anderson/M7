//
// Created by rja on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"

/**
 * specification of the extents and properties of the fermionic single particle basis
 */
struct FrmBasisData {
    /**
     * number of electrons in the system
     */
    const size_t m_nelec;
    /**
     * number of sites or orbitals
     */
    const size_t m_nsite;
    /**
     * number of fermionic degrees of freedom (2*nsite)
     */
    const size_t m_nspinorb;
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

    FrmBasisData(size_t nelec, size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved, int ms2, bool ms2_conserve);
    /*
     * default ctor: non-resolved spin, and no Ms2 restriction
     */
    FrmBasisData(size_t nelec, size_t nsite):
        FrmBasisData(nelec, nsite, {nsite}, false, 0, false){}

    bool operator==(const FrmBasisData& other) const;

    operator bool() const {
        return m_nsite;
    }
};

/**
 * specification of the extents and properties of the bosonic single particle basis
 */
struct BosBasisData {
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

    BosBasisData(size_t nmode, size_t nboson, bool nboson_conserve, size_t occ_cutoff):
            m_nmode(nmode), m_nboson(nboson), m_nboson_conserve(nboson_conserve), m_occ_cutoff(occ_cutoff){}
    BosBasisData(size_t nmode): BosBasisData(nmode, 0ul, false, defs::max_bos_occ){}
    bool operator==(const BosBasisData& other) const;

    operator bool() const {
        return m_nmode;
    }
};

/**
 * admits a common interface for all initialisations of MBFs
 */
struct BasisData {
    const FrmBasisData m_frm;
    const BosBasisData m_bos;

    BasisData(FrmBasisData frm, BosBasisData bos);
    BasisData(size_t nsite, size_t nmode): m_frm(nsite), m_bos(nmode){}
    BasisData(): BasisData(0,0){}

    void require_pure_frm() const;
    void require_pure_bos() const;
    bool operator==(const BasisData& other) const;

    operator bool() const {
        return m_frm || m_bos;
    }
};

#endif //M7_BASISDATA_H