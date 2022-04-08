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
     * mapping from fermion site indices to Abelian group labels
     */
    const AbelianGroupMap m_abgrp_map;

    FrmBasisData(size_t nsite, AbelianGroupMap abgrp_map);
    explicit FrmBasisData(size_t nsite);
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
    /**
     * maximum allowed occupation of each bosonic mode
     */
    const defs::bos_occ_t m_nboson_max;

    BosBasisData(size_t nmode, size_t nboson_max);
    explicit BosBasisData(size_t nmode): BosBasisData(nmode, defs::max_bos_occ){}
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

    void require_pure_frm() const;
    void require_pure_bos() const;
    bool operator==(const BasisData& other) const;
};


#endif //M7_BASISDATA_H
