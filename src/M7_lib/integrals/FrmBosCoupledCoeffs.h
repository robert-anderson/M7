//
// Created by Robert J. Anderson on 28/08/2021.
//

#ifndef M7_FRMBOSCOUPLEDCOEFFS_H
#define M7_FRMBOSCOUPLEDCOEFFS_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/parallel/SharedArray.h>

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
    uint_t index(uint_t n, uint_t p, uint_t q) const;

public:
    /**
     * size of the many-body bases
     */
    const sys::Size m_sizes;
    /**
     * the extent of the fermion indices in the array and its square to give the stride between consecutive mode indices
     */
    const uint_t m_ncoeff_ind_frm, m_ncoeff_ind_frm2;
    /**
     * the extent of the boson indices (always the number of modes since boson indices are never spin resolved)
     */
    const uint_t m_ncoeff_ind_bos;
    /**
     * array of "integrals" only stored on the root-rank of each node
     */
    SharedArray<defs::ham_t> m_v;

    FrmBosCoupledCoeffs(sys::Size sizes, bool spin_resolved);

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
    void set(uint_t n, uint_t p, uint_t q, defs::ham_t value);

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
    defs::ham_t get(uint_t n, uint_t p, uint_t q) const;

};


#endif //M7_FRMBOSCOUPLEDCOEFFS_H
