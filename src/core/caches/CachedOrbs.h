//
// Created by rja on 24/08/2021.
//

#ifndef M7_CACHEDORBS_H
#define M7_CACHEDORBS_H

#include <src/core/caches/DecodedDeterminants.h>

struct FrmOnvField;
struct BosOnvField;
struct FrmBosOnvField;

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbital indices when
 * needed
 */
namespace decoded_mbf {

    namespace spinorbs {

        /**
         * base class for the unstructured decoding of a FrmOnv (bit representation) into a simple array of set or clear 
         * bit indices.
         * the subclasses must provide the logic for updating based on the clear or set positions.
         */
        struct SimpleBase {
            /**
             * spin orbital indices
             */
            defs::inds m_inds;

        public:
            size_t size() const;

            const size_t &operator[](const size_t &i) const;

            const defs::inds &inds() const;

            void clear();

            bool empty();
        };

        /**
         * update based on the set positions
         */
        struct SimpleOccs : SimpleBase {
            void update(const FrmOnvField &mbf);
        };

        /**
         * update based on the clr positions
         */
        struct SimpleVacs : SimpleBase {
            void update(const FrmOnvField &mbf);
        };


        /**
         * base class for the structured decoding of a FrmOnv (bit representation) into a ragged array of set or clear
         * bit indices segregated by the label assigned to the spinorbital via the given integer map
         * the subclasses must provide the logic for updating based on the clear or set positions.
         */
        struct LabelledBase {
            /**
             * ragged array of vectors to store set or clear positions. one vector per label
             */
            std::vector<defs::inds> m_inds;
            /**
             * map from the spinorb index to the "label"
             */
            const defs::inds m_map;

            /**
             * duplicate the spatial orbital irrep label map into two spin channels (spin major mapping)
             * @param grp_map
             *  spatial orbital irrep map
             * @return
             *  spin orbital irrep map with 2*nirrep labels
             */
            static defs::inds make_spinorb_map(const defs::inds &site_irreps, size_t nirrep);

            LabelledBase(size_t nelement, const defs::inds &map);

            LabelledBase(const LabelledBase &other);

            LabelledBase(LabelledBase &&other);

            LabelledBase &operator=(const LabelledBase &other);

            LabelledBase &operator=(LabelledBase &&other);

            size_t size(const size_t &ielement) const;

            /**
             * @param i
             *  label index
             * @return
             *  reference to the set or clear indices for label index i
             */
            const defs::inds &operator[](const size_t &i) const;

            void clear();
        };

        /**
         * update based on the set positions
         */
        struct LabelledOccs : LabelledBase {
            /**
             * flat orbitals are included so they can be decoded in the same loop as the labelled indices
             */
            SimpleOccs m_simple;
        protected:
            LabelledOccs(size_t nelement, const defs::inds &map);

        public:
            void update(const FrmOnvField &mbf);

            void clear();

            bool empty();
        };

        /**
         * update based on the clr positions
         */
        struct LabelledVacs : LabelledBase {
            /**
             * flat orbitals are included so they can be decoded in the same loop as the labelled indices
             */
            SimpleVacs m_simple;

        protected:
            LabelledVacs(size_t nelement, const defs::inds &map);

        public:
            void update(const FrmOnvField &mbf);

            void clear();

            bool empty();
        };


        /**
         * multidimensionality is introduced in such a way that the decoding of the above classes is agnostic to the
         * shape of the label array, and the multidimensional structure defined here is agnostic to the details of
         * determinantal decoding, since the flat index array is only referred to by a const ref.
         * @tparam nind
         *  number of elements in shape of label array
         */
        template<size_t nind>
        struct NdBase {
            const NdFormat<nind> m_format;
        protected:
            /**
             * ragged array of indices from a LabelledOccs or LabelledVacs instance
             */
            const std::vector<defs::inds> &m_inds_ref;

        public:
            NdBase(std::array<size_t, nind> shape, const std::vector<defs::inds> &inds) :
                    m_format(shape), m_inds_ref(inds) {}

            size_t size(const std::array<size_t, nind> &inds) const {
                return m_inds_ref[m_format.flatten(inds)].size();
            }

            const defs::inds &operator[](const std::array<size_t, nind> &inds) const {
                return m_inds_ref[m_format.flatten(inds)];
            }
        };

        /**
         * label sets will often make sense as multidimensional objects, here the flat decoding of determinantal bit
         * representations (LabelledOccs) is combined with the multidimensionality (NdBase)
         * @tparam nind
         *  number of elements in shape of label array
         */
        template<size_t nind>
        struct NdLabelledOccs : LabelledOccs, NdBase<nind> {
            using LabelledOccs::size;
            using NdBase<nind>::size;
            using LabelledOccs::operator[];
            using NdBase<nind>::operator[];

            NdLabelledOccs(std::array<size_t, nind> shape, const defs::inds &map) :
                    LabelledOccs(NdFormat<nind>(shape).m_nelement, map), NdBase<nind>(shape, m_inds) {}
        };

        /**
         * label sets will often make sense as multidimensional objects, here the flat decoding of determinantal bit
         * representations (LabelledVacs) is combined with the multidimensionality (NdBase)
         * @tparam nind
         *  number of elements in shape of label array
         */
        template<size_t nind>
        struct NdLabelledVacs : LabelledVacs, NdBase<nind> {
            using LabelledVacs::size;
            using NdBase<nind>::size;
            using LabelledVacs::operator[];
            using NdBase<nind>::operator[];

            NdLabelledVacs(std::array<size_t, nind> shape, const defs::inds &map) :
                    LabelledVacs(NdFormat<nind>(shape).m_nelement, map), NdBase<nind>(shape, m_inds) {}
        };

        /**
         * The immediate practical application of the general multidimensional decoding defined above, is in the case
         * where there is utility in segregating the set or clear indices by both spin and point group irrep.
         */
        struct SpinSymOccs : NdLabelledOccs<2> {
            explicit SpinSymOccs(const AbelianGroupMap &grp_map) :
                    NdLabelledOccs<2>({2, grp_map.m_grp.nirrep()},
                                      make_spinorb_map(grp_map.m_site_irreps, grp_map.m_grp.nirrep())) {}
        };

        /**
         * vacant orbital analogue of the above definition
         */
        struct SpinSymVacs : NdLabelledVacs<2> {
            explicit SpinSymVacs(const AbelianGroupMap &grp_map) :
                    NdLabelledVacs<2>({2, grp_map.m_grp.nirrep()},
                                      make_spinorb_map(grp_map.m_site_irreps, grp_map.m_grp.nirrep())) {}
        };
    }

    class FrmOnv {
        const FrmOnvField &m_mbf;
        /**
         * spin/sym-partitioned occupied spin orbitals
         */
        spinorbs::SpinSymOccs m_occ;
        /**
         * spin/sym-partitioned vacant spin orbitals
         */
        spinorbs::SpinSymVacs m_vac;
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
        FrmOnv(const FrmOnvField &mbf, const AbelianGroupMap &grp_map = {});

        /**
         * clear all cached assets
         */
        void clear();

        /**
         * call update on the occupied orbitals if required
         * @return
         *  the occupied spin orbitals
         */
        const spinorbs::SpinSymOccs &occ();

        /**
         * call update on the occupied orbitals if required
         * @return
         *  the occupied spin orbitals
         */
        const spinorbs::SpinSymVacs &vac();

        /**
         * update the vector of labels of nonempty occupied and vacant spin orbs
         * @return
         *  vector of labels of nonempty occupied and vacant spin orbs
         */
        const defs::inds &nonempty_pair_labels();

        /**
         * update the vector of site indices with at least one fermion
         * @return
         *  vector of site indices with at least one fermion
         */
        const defs::inds &occ_sites();

        /**
         * update the vector of site indices occupied with two fermions
         * @return
         *  vector of site indices with two fermions
         */
        const defs::inds &doubly_occ_sites();

        /**
         * update the vector of site indices occupied with zero or two fermions
         * @return
         *  vector of site indices with zero or two fermions
         */
        const defs::inds &not_singly_occ_sites();
    };

    class BosOnv {
        const BosOnvField &m_mbf;
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

    public:
        BosOnv(const BosOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();

        /**
         * update the vector mode indices of boson creation operators which would create mbf upon application to the vacuum
         * @param mbf
         *  bosonic many-body basis function
         * @return
         *  vector of occupied bosonic modes with repetition
         */
        const defs::inds &bos_op_inds();

        const defs::inds &occ_bos_inds();
    };

    class FrmBosOnv {
        const FrmBosOnvField &m_mbf;
        /**
         * site indices with any electrons and any associated bosons (imode==isite) e.g. holstein
         */
        defs::inds m_occ_sites_nonzero_bosons;

    public:
        FrmBosOnv(const FrmBosOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();

        /**
         * update the vector of site indices with at least one fermion and at least one boson
         * @param mbf
         *  many-body basis function
         * @return
         *  vector of site indices with at least one fermion and at least one boson
         */
        const defs::inds &occ_sites_nonzero_bosons();
    };
};


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
    CachedOrbs(const AbelianGroupMap &grp_map = {});

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

#endif //M7_CACHEDORBS_H