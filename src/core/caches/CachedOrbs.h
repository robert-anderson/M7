//
// Created by rja on 24/08/2021.
//

#ifndef M7_CACHEDORBS_H
#define M7_CACHEDORBS_H

#include <algorithm>
#include <src/core/basis/AbelianGroup.h>
#include "src/core/parallel/MPIAssert.h"
#include "src/core/nd/NdFormat.h"
#include "src/defs.h"

struct FrmOnvField;
struct BosOnvField;
struct FrmBosOnvField;

/**
 * provides a store which only updates the flat and spin/sym-partitioned occupied and vacant spin orbital indices when
 * needed
 */
namespace decoded_mbf {

    /**
     * base class for MBF "caches", i.e. simple or structured arrays updated based on the current value of an associated
     * many body basis function object
     */
    struct Cache {
        /**
         * only used for debugging as a check that the current state of the cache corresponds to that of the MBF object
         * @return
         *  true if the current hash matches that of the last update
         */
        virtual bool is_valid() const = 0;
    protected:
        /**
         * should be updated only in the debug build any time an update is performed
         */
        defs::hash_t m_last_update_hash = 0;
    };

    /**
     * in most cases, a cache only requires a flat integer vector in which to store indices
     */
    struct SimpleContainer {
    protected:
        defs::inds m_inds;

    public:

        void clear();

        bool empty();
    };

    namespace frm {

        struct Base : Cache {
        protected:
            const FrmOnvField &m_mbf;
        public:
            Base(const FrmOnvField &mbf);

            bool is_valid() const override;
        };

        /**
         * base class for the unstructured decoding of a FrmOnv (bit representation) into a simple array of set or clear 
         * bit indices.
         * the subclasses must provide the logic for updating based on the clear or set positions.
         */
        struct SimpleBase : Base, SimpleContainer {
        protected:
            const defs::inds& validated() const;
        public:
            explicit SimpleBase(const FrmOnvField &mbf);
        };

        /**
         * update based on the set positions
         */
        struct SimpleOccs : SimpleBase {
            explicit SimpleOccs(const FrmOnvField &mbf);
            const defs::inds& get();
        };

        /**
         * update based on the clr positions
         */
        struct SimpleVacs : SimpleBase {
            explicit SimpleVacs(const FrmOnvField &mbf);
            const defs::inds& get();
        };


        /**
         * base class for the structured decoding of a FrmOnv (bit representation) into a ragged array of set or clear
         * bit indices segregated by the label assigned to the spinorbital via the given integer map
         * the subclasses must provide the logic for updating based on the clear or set positions.
         */
        struct LabelledBase : Base {
        protected:
            /**
             * ragged array of vectors to store set or clear positions. one vector per label
             */
            std::vector<defs::inds> m_inds;
            /**
             * map from the spinorb index to the "label"
             */
            const defs::inds m_map;
            /**
             * the simple decoding is cached at the same time, since it is often useful and available at small
             * additional cost
             */
             defs::inds m_simple_inds;

            LabelledBase(size_t nelement, const defs::inds &map, const FrmOnvField &mbf);

            const std::vector<defs::inds>& validated() const;

        public:
            /**
             * duplicate the spatial orbital irrep label map into two spin channels (spin major mapping)
             * @param grp_map
             *  spatial orbital irrep map
             * @return
             *  spin orbital irrep map with 2*nirrep labels
             */
            static defs::inds make_spinorb_map(const defs::inds &site_irreps, size_t nirrep);

            void clear();

            bool empty();
        };

        /**
         * update based on the set positions
         */
        struct LabelledOccs : LabelledBase {
        protected:
            LabelledOccs(size_t nelement, const defs::inds &map, const FrmOnvField &mbf);

        public:
            const std::vector<defs::inds>& get();
            const defs::inds& simple();
        };

        /**
         * update based on the clr positions
         */
        struct LabelledVacs : LabelledBase {
        protected:
            LabelledVacs(size_t nelement, const defs::inds &map, const FrmOnvField &mbf);

        public:
            const std::vector<defs::inds>& get();
            const defs::inds& simple();
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

            size_t size(const size_t &i) const {
                return m_inds_ref[i].size();
            }
            size_t size(const std::array<size_t, nind> &inds) const {
                return m_inds_ref[m_format.flatten(inds)].size();
            }
            const defs::inds &operator[](const size_t &i) const {
                return m_inds_ref[i];
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
        struct NdLabelledOccs : LabelledOccs {
        protected:
            NdBase<nind> m_nd_inds;
        public:
            NdLabelledOccs(std::array<size_t, nind> shape, const defs::inds &map, const FrmOnvField &mbf) :
                    LabelledOccs(NdFormat<nind>(shape).m_nelement, map, mbf),
                    m_nd_inds(shape, m_inds) {}
            const NdBase<nind>& get(){
                LabelledOccs::get();
                return m_nd_inds;
            }

        protected:
            const NdBase<nind>& validated() const {
                DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
                return m_nd_inds;
            }
        };

        /**
         * label sets will often make sense as multidimensional objects, here the flat decoding of determinantal bit
         * representations (LabelledVacs) is combined with the multidimensionality (NdBase)
         * @tparam nind
         *  number of elements in shape of label array
         */
        template<size_t nind>
        struct NdLabelledVacs : LabelledVacs {
        protected:
            NdBase<nind> m_nd_inds;
        public:
            NdLabelledVacs(std::array<size_t, nind> shape, const defs::inds &map, const FrmOnvField &mbf) :
                    LabelledVacs(NdFormat<nind>(shape).m_nelement, map, mbf),
                    m_nd_inds(shape, m_inds) {}
            const NdBase<nind>& get(){
                LabelledVacs::get();
                return m_nd_inds;
            }

        protected:
            const NdBase<nind>& validated() const {
                DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
                return m_nd_inds;
            }
        };

        /**
         * The immediate practical application of the general multidimensional decoding defined above, is in the case
         * where there is utility in segregating the set or clear indices by both spin and point group irrep.
         */
        struct SpinSymOccs : NdLabelledOccs<2> {
            SpinSymOccs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf);
        };

        /**
         * vacant orbital analogue of the above definition
         */
        struct SpinSymVacs : NdLabelledVacs<2> {
            SpinSymVacs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf);
        };

        /**
         * labels (flat indices of the spin-sym partitioning) with at least one occupied and one vacant orbital
         */
        struct NonEmptyPairLabels : SimpleBase {
        protected:
            SpinSymOccs &m_occs;
            SpinSymVacs &m_vacs;
        public:
            NonEmptyPairLabels(const FrmOnvField &mbf, SpinSymOccs &occs, SpinSymVacs &vacs);

            const defs::inds& get();
        };

    }

    struct FrmOnv {
        /**
         * occupied spin orbitals
         */
        frm::SimpleOccs m_simple_occs;
        /**
         * vacant spin orbitals
         */
        frm::SimpleVacs m_simple_vacs;
        /**
         * spin/sym-partitioned occupied spin orbitals
         */
        frm::SpinSymOccs m_spin_sym_occs;
        /**
         * spin/sym-partitioned vacant spin orbitals
         */
        frm::SpinSymVacs m_spin_sym_vacs;
        /**
         * labels (flat indices of the spin-sym partitioning) with at least one occupied and one vacant orbital
         */
        frm::NonEmptyPairLabels m_nonempty_pair_labels;
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

        FrmOnv(const FrmOnvField &mbf, const AbelianGroupMap &grp_map);

        FrmOnv(const FrmOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();

#if 0
        /**
         * call update on the simply-decoded occupied orbitals if required
         * @return
         *  the occupied spin orbitals
         */
        const frm::SimpleOccs &simple_occs();

        /**
         * call update on the simply-decoded vacant orbitals if required
         * @return
         *  the vacant spin orbitals
         */
        const frm::SimpleVacs &simple_vacs();

        /**
         * call update on the spin and point group symmetry-segregated occupied orbitals if required
         * @return
         *  the occupied spin orbitals
         */
        const frm::SpinSymOccs &spin_sym_occs();

        /**
         * call update on the spin and point group symmetry-segregated vacant orbitals if required
         * @return
         *  the vacant spin orbitals
         */
        const frm::SpinSymVacs &spin_sym_vacs();

        /**
         * update the vector of labels of nonempty occupied and vacant spin orbs
         * @return
         *  vector of labels of nonempty occupied and vacant spin orbs
         */
        const defs::inds &nonempty_pair_labels();

#endif
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

#endif //M7_CACHEDORBS_H