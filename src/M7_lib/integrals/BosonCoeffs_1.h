//
// Created by Robert J. Anderson on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_1_H
#define M7_BOSONCOEFFS_1_H

#include <M7_lib/parallel/SharedArray.h>

class BosonCoeffs_1 {

    uint_t index(uint_t n, uint_t m) const {
        return n * m_nmode + m;
    }

public:
    const uint_t m_nmode;
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_1(uint_t nmode);

    void set(uint_t n, uint_t m, defs::ham_t value);

    defs::ham_t get(uint_t n, uint_t m) const;
};


#endif //M7_BOSONCOEFFS_1_H
