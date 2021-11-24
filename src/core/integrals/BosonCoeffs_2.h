//
// Created by rja on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_2_H
#define M7_BOSONCOEFFS_2_H

#include "src/core/parallel/SharedArray.h"
#include "Integrals.h"

class BosonCoeffs_2 {

    size_t index(const size_t &i, const size_t &j, const size_t& k, const size_t& l) const {
        auto ij = i<=j ? trig(i, j) : trig(j, i);
        auto kl = k<=l ? trig(k, l) : trig(l, k);
        return ij<=kl ? trig(ij, kl) : trig(kl, ij);
    }

public:
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_2(size_t nmode);

    void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;

    const defs::ham_t& phys_element(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;
};


#endif //M7_BOSONCOEFFS_2_H
