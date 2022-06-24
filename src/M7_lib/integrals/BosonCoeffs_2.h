//
// Created by Robert J. Anderson on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_2_H
#define M7_BOSONCOEFFS_2_H

#include <M7_lib/parallel/SharedArray.h>
#include <M7_lib/util/Integer.h>

class BosonCoeffs_2 {

    uint_t index(uint_t i, uint_t j, uint_t k, uint_t l) const {
        using namespace integer;
        auto ij = i*m_nmode+j;
        auto kl = k*m_nmode+l;
        return ij>=kl ? trigmap(ij, kl) : trigmap(kl, ij);
    }

    const uint_t m_nmode;
public:
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_2(uint_t nmode);

    void set(uint_t i, uint_t j, uint_t k, uint_t l, defs::ham_t value);

    defs::ham_t get(uint_t i, uint_t j, uint_t k, uint_t l) const;

    defs::ham_t phys_element(uint_t i, uint_t j, uint_t k, uint_t l) const;
};


#endif //M7_BOSONCOEFFS_2_H
