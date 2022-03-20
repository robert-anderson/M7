//
// Created by rja on 14/09/2021.
//

#ifndef M7_BASISDATA_H
#define M7_BASISDATA_H

#include <utility>

#include <M7_lib/parallel/MPIAssert.h>
#include "AbelianGroup.h"

/**
 * admits a common interface for all initialisations of MBFs
 */
struct BasisData {
    /**
     * number of fermionic sites or orbitals
     */
    size_t m_nsite;
    /**
     * number of bosonic modes
     */
    size_t m_nmode;
    /**
     * number of fermionic degrees of freedom (2*nsite)
     */
    size_t m_nspinorb;
    /**
     * mapping from fermion site indices to abelian group labels
     */
    AbelianGroupMap m_frm_abgrp_map;
    BasisData(size_t nsite, size_t nmode, AbelianGroupMap frm_abgrp_map):
        m_nsite(nsite), m_nmode(nmode), m_nspinorb(nsite*2), m_frm_abgrp_map(std::move(frm_abgrp_map)){}

    BasisData(size_t nsite, size_t nmode): BasisData(nsite, nmode, {nsite}){}

    void require_pure_frm() const {
        REQUIRE_FALSE(m_nmode, "MBF specification is not purely fermionic");
    }
    void require_pure_bos() const {
        REQUIRE_FALSE(m_nsite, "MBF specification is not purely bosonic");
    }
};


#endif //M7_BASISDATA_H
