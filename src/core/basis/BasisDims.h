//
// Created by rja on 14/09/2021.
//

#ifndef M7_BASISDIMS_H
#define M7_BASISDIMS_H

#include "src/core/parallel/MPIAssert.h"
/**
 * admits a common interface for all initialisations of MBFs
 */
struct BasisDims {
    size_t m_nsite;
    size_t m_nmode;
    size_t m_nspinorb;
    BasisDims(size_t nsite, size_t nmode):
        m_nsite(nsite), m_nmode(nmode), m_nspinorb(nsite*2){}

    operator bool() const {
        return m_nsite || m_nmode;
    }

    void require_pure_frm() {
        REQUIRE_FALSE(m_nmode, "MBF specification is not purely fermionic");
    }
    void require_pure_bos() {
        REQUIRE_FALSE(m_nsite, "MBF specification is not purely bosonic");
    }
};


#endif //M7_BASISDIMS_H
