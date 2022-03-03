//
// Created by rja on 14/09/2021.
//

#ifndef M7_BASISDIMS_H
#define M7_BASISDIMS_H

#include <utility>

#include "src/core/parallel/MPIAssert.h"
#include "AbelianGroup.h"

/**
 * admits a common interface for all initialisations of MBFs
 */
struct BasisDims {
    size_t m_nsite;
    size_t m_nmode;
    size_t m_nspinorb;
    AbelianGroupMap m_frm_abgrp_map;
    BasisDims(size_t nsite, size_t nmode, AbelianGroupMap frm_abgrp_map):
        m_nsite(nsite), m_nmode(nmode), m_nspinorb(nsite*2), m_frm_abgrp_map(std::move(frm_abgrp_map)){}

    BasisDims(size_t nsite, size_t nmode): BasisDims(nsite, nmode, {nsite}){}

    void require_pure_frm() const {
        REQUIRE_FALSE(m_nmode, "MBF specification is not purely fermionic");
    }
    void require_pure_bos() const {
        REQUIRE_FALSE(m_nsite, "MBF specification is not purely bosonic");
    }
};


#endif //M7_BASISDIMS_H
