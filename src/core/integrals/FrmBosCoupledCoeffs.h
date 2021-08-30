//
// Created by rja on 28/08/2021.
//

#ifndef M7_FRMBOSCOUPLEDCOEFFS_H
#define M7_FRMBOSCOUPLEDCOEFFS_H

#include "src/core/parallel/SharedArray.h"

class FrmBosCoupledCoeffs {

    size_t index(const size_t &n, const size_t &p, const size_t &q) const;

public:
    const size_t m_nmode, m_nmode2;
    SharedArray<defs::ham_t> m_v;

    FrmBosCoupledCoeffs(size_t nmode);

    void set(const size_t& n, const size_t& p, const size_t& q, const defs::ham_t& value);

    const defs::ham_t& get(const size_t& n, const size_t& p, const size_t& q) const;

    bool constant_diagonal() const;
};


#endif //M7_FRMBOSCOUPLEDCOEFFS_H
