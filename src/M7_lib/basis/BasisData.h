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

    FrmBasisData(size_t nsite, AbelianGroupMap abgrp_map, bool spin_resolved);
    FrmBasisData(size_t nsite);
    FrmBasisData();
    bool operator==(const FrmBasisData& other) const;
};

/**
 * specification of the extents and properties of the bosonic single particle basis
 */
struct BosBasisData {
    /**
     * number of bosonic modes
     */
    const size_t m_nmode;

    BosBasisData(size_t nmode);
    BosBasisData(): BosBasisData(0ul){}
    bool operator==(const BosBasisData& other) const;
};

/**
 * admits a common interface for all initialisations of MBFs
 */
struct BasisData {
    const FrmBasisData m_frm;
    const BosBasisData m_bos;

    BasisData(FrmBasisData frm, BosBasisData bos);
    BasisData(size_t nsite, size_t nmode): m_frm(nsite), m_bos(nmode){}

    void require_pure_frm() const;
    void require_pure_bos() const;
    bool operator==(const BasisData& other) const;
};

#endif //M7_BASISDATA_H