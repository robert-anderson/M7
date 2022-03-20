//
// Created by rja on 02/03/2022.
//

#ifndef M7_CACHEDORBSOLD_H
#define M7_CACHEDORBSOLD_H

#include "DecodedDeterminants.h"

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbital indices when
 * needed
 */
class CachedOrbs {
    /**
     * spin/sym-partitioned occupied spin orbitals
     */
    SpinSymOccOrbs__ m_occ;
    /**
     * spin/sym-partitioned vacant spin orbitals
     */
    SpinSymVacOrbs__ m_vac;
    /**
     * labels (flat indices of the partitioning) with at least one occupied and one vacant orbital
     */
    defs::inds m_nonempty_pair_labels;
    /**
     * site indices with any electrons
     */
    defs::inds m_occ_sites;
    /**
     * site indices with any electrons and any associated bosons (imode==isite) e.g. holstein
     */
    defs::inds m_occ_sites_nonzero_bosons;
    /**
     * site indices with 2 electrons
     */
    defs::inds m_doubly_occ_sites;
    /**
     * site indices with 0 or 2 electrons
     */
    defs::inds m_not_singly_occ_sites;
    /**
     * boson mode indices with repetition
     *  e.g. [0, 2, 0, 3, 1] decodes as:
     *      [1, 1, 3, 3, 3, 4]
     */
    defs::inds m_bos_op_inds;
    /**
     * indices of boson modes with any non-zero occupation
    *  e.g. [0, 2, 0, 3, 1] decodes as:
    *      [1, 3, 4]
     */
    defs::inds m_occ_bos_inds;

    /**
     * clear only those cached assets which depend on both the fermion and boson sectors of the src MBF
     */
    void clear_frmbos_only();

    /**
     * clear all fermion sector-dependent cached assets
     */
    void clear_frm();

    /**
     * clear all boson sector-dependent cached assets
     */
    void clear_bos();

public:
    CachedOrbs(const AbelianGroupMap &grp_map);

    CachedOrbs(size_t nsite);

    void clear(const field::FrmOnv &mbf);

    void clear(const field::BosOnv &mbf);

    void clear(const field::FrmBosOnv &mbf);

    /**
     * alias for the FrmBosOnv overload: clear all cached assets
     */
    void clear();

    /**
     * call update on the occupied orbitals if required
     * @param mbf
     *  many-body basis function
     * @return
     *  the occupied spin orbitals
     */
    const SpinSymOccOrbs__ &occ(const field::FrmOnv &mbf);

    /**
     * call update on the occupied orbitals if required
     * @param mbf
     *  many-body basis function
     * @return
     *  the occupied spin orbitals
     */
    const SpinSymVacOrbs__ &vac(const field::FrmOnv &mbf);

    /**
     * update the vector of labels of nonempty occupied and vacant spin orbs
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of labels of nonempty occupied and vacant spin orbs
     */
    const defs::inds &nonempty_pair_labels(const field::FrmOnv &mbf);

    /**
     * update the vector of site indices with at least one fermion
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with at least one fermion
     */
    const defs::inds &occ_sites(const field::FrmOnv &mbf);

    /**
     * update the vector of site indices with at least one fermion and at least one boson
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with at least one fermion and at least one boson
     */
    const defs::inds &occ_sites_nonzero_bosons(const field::FrmBosOnv &mbf);

    /**
     * update the vector of site indices occupied with two fermions
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with two fermions
     */
    const defs::inds &doubly_occ_sites(const field::FrmOnv &mbf);

    /**
     * update the vector of site indices occupied with zero or two fermions
     * @param mbf
     *  many-body basis function
     * @return
     *  vector of site indices with zero or two fermions
     */
    const defs::inds &not_singly_occ_sites(const field::FrmOnv &mbf);

    /**
     * update the vector mode indices of boson creation operators which would create mbf upon application to the vacuum
     * @param mbf
     *  bosonic many-body basis function
     * @return
     *  vector of occupied bosonic modes with repetition
     */
    const defs::inds &bos_op_inds(const field::BosOnv &mbf);

    const defs::inds &occ_bos_inds(const field::BosOnv &mbf);
};


#endif //M7_CACHEDORBSOLD_H
