//
// Created by rja on 02/09/2021.
//

#ifndef M7_BOSONCOEFFS_H
#define M7_BOSONCOEFFS_H

#include "src/core/parallel/SharedArray.h"

class BosonCoeffs {

    size_t index(const size_t &n, const size_t &m) const {
        return n * m_nmode + m;
    }

public:
    const size_t m_nmode;
    SharedArray<defs::ham_t> m_v;

    BosonCoeffs(size_t nmode);

    void set(const size_t& n, const size_t& m, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& n, const size_t& m) const;
};


#endif //M7_BOSONCOEFFS_H
