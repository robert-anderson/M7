//
// Created by rja on 28/08/2021.
//

#ifndef M7_FRMBOSCOUPLEDCOEFFS_H
#define M7_FRMBOSCOUPLEDCOEFFS_H

#include <src/core/basis/BasisData.h>
#include "src/core/parallel/SharedArray.h"

class FrmBosCoupledCoeffs {

    /**
     * @param n
     *  mode index of boson creation operator
     * @param p
     *  spin orbital index of fermion creation operator
     * @param q
     *  spin orbital index of fermion annihilation operator
     * @return
     *  flat index of the element in the array
     */
    size_t index(size_t n, size_t p, size_t q) const;

public:
    const BasisData m_bd;
    /**
     * the extent of the fermion indices in the array and its square to give the stride between consecutive mode indices
     */
    const size_t m_nintind_frm, m_nintind_frm2;
    /**
     * the extent of the boson indices (always the number of modes since boson indices are never spin resolved)
     */
    const size_t m_nintind_bos;
    /**
     * array of "integrals" only stored on the root-rank of each node
     */
    SharedArray<defs::ham_t> m_v;

    FrmBosCoupledCoeffs(BasisData bd, bool spin_resolved);

    /**
     * assign a value to the indexed element (should only be called on the root rank of each node)
     * @param n
     *  mode index of boson creation operator
     * @param p
     *  spin orbital index of fermion creation operator
     * @param q
     *  spin orbital index of fermion annihilation operator
     * @param value
     *  coefficient value
     */
    void set(const size_t& n, const size_t& p, const size_t& q, const defs::ham_t& value);

    /**
     * access the coefficient value (called on any rank)
     * @param n
     *  mode index of boson creation operator
     * @param p
     *  spin orbital index of fermion creation operator
     * @param q
     *  spin orbital index of fermion annihilation operator
     * @return
     *  coefficient value
     */
    const defs::ham_t& get(const size_t& n, const size_t& p, const size_t& q) const;
    /**
     * only applies in cases where the number of sites is equal to number of modes
     * @return
     *  true if applicable and all V[i, i, i] have the same value (e.g. Hubbard--Holstein)
     */
    bool constant_diagonal() const;
};


#endif //M7_FRMBOSCOUPLEDCOEFFS_H
