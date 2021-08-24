//
// Created by rja on 24/08/2021.
//

#ifndef M7_CACHEDORBS_H
#define M7_CACHEDORBS_H

#include <src/core/basis/DecodedDeterminant.h>

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbtial indices when
 * needed
 */
class CachedOrbs {
    /**
     * spin/sym-partitioned occupied spin orbitals
     */
    SpinSymOccOrbs m_occ;
    /**
     * spin/sym-partitioned vacant spin orbitals
     */
    SpinSymVacOrbs m_vac;
    /**
     * labels (flat indices of the partitioning) with at least one occupied and one vacant orbital
     */
    defs::inds m_nonempty_pair_labels;
public:
    CachedOrbs(const AbelianGroupMap& grp_map);
    /**
     * set both occ and vac spin orbital stores to "un-updated" status
     */
    void clear();
    /**
     * call update on the occupied orbitals if required
     * @param mbf
     *  many-body basis function
     * @return
     *  the occupied spin orbitals
     */
    const SpinSymOccOrbs& occ(const field::FrmOnv &mbf);
    /**
     * call update on the occupied orbitals if required
     * @param mbf
     *  many-body basis function
     * @return
     *  the occupied spin orbitals
     */
    const SpinSymVacOrbs& vac(const field::FrmOnv &mbf);
    /**
     * update the vector of labels of nonempty occupied and vacant spin orbs
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of labels of nonempty occupied and vacant spin orbs
     */
    const defs::inds& nonempty_pair_labels(const field::FrmOnv &mbf);
};


#endif //M7_CACHEDORBS_H
