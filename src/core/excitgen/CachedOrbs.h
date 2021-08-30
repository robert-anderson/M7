//
// Created by rja on 24/08/2021.
//

#ifndef M7_CACHEDORBS_H
#define M7_CACHEDORBS_H

#include <src/core/basis/DecodedDeterminants.h>

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbital indices when
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
    /**
     * site indices with any electrons
     */
    defs::inds m_occ_sites;
    /**
     * site indices with 2 electrons
     */
    defs::inds m_doubly_occ_sites;
    /**
     * site indices with 0 or 2 electrons
     */
    defs::inds m_not_singly_occ_sites;
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
    /**
     * update the vector of site indices with at least one fermion
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with at least one fermion
     */
    const defs::inds& occ_sites(const field::FrmOnv &mbf);
    /**
     * update the vector of site indices occupied with two fermions
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with two fermions
     */
    const defs::inds& doubly_occ_sites(const field::FrmOnv &mbf);
    /**
     * update the vector of site indices occupied with zero or two fermions
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with zero or two fermions
     */
    const defs::inds& not_singly_occ_sites(const field::FrmOnv &mbf);
};


#endif //M7_CACHEDORBS_H
