//
// Created by Robert J. Anderson on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_1_H
#define M7_BOSONCOEFFS_1_H

#include <M7_lib/parallel/SharedArray.h>

class BosonCoeffs_1 {

    size_t index(const size_t &n, const size_t &m) const {
        return n * m_nmode + m;
    }

public:
    const size_t m_nmode;
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs_1(size_t nmode);

    void set(const size_t& n, const size_t& m, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& n, const size_t& m) const;
};


#endif //M7_BOSONCOEFFS_1_H
