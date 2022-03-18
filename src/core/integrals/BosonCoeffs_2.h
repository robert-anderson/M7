//
// Created by rja on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_2_H
#define M7_BOSONCOEFFS_2_H

#include "parallel/SharedArray.h"

#include "Integrals.h"

class BosonCoeffs_2 {

    size_t index(const size_t &i, const size_t &j, const size_t& k, const size_t& l) const {
        auto ij = i*m_nmode+j;
        auto kl = k*m_nmode+l;
        return ij<=kl ? trig(ij, kl) : trig(kl, ij);
    }

    const size_t m_nmode;
public:
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_2(size_t nmode);

    void set(const size_t& i, const size_t& j, const size_t& k, const size_t& l, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;

    const defs::ham_t& phys_element(const size_t& i, const size_t& j, const size_t& k, const size_t& l) const;
};


#endif //M7_BOSONCOEFFS_2_H
