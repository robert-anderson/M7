//
// Created by Robert J. Anderson on 04/03/2022.
//

#ifndef M7_DECODEDFRMONV_H
#define M7_DECODEDFRMONV_H

#include "DecodedMbf.h"

namespace decoded_mbf {
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
            const defs::inds_t &validated() const;

        public:
            explicit SimpleBase(const FrmOnvField &mbf);
        };

        /**
         * update based on the set positions
         */
        struct SimpleOccs : SimpleBase {
            explicit SimpleOccs(const FrmOnvField &mbf);

            const defs::inds_t &get();
        };

        /**
         * update based on the clr positions
         */
        struct SimpleVacs : SimpleBase {
            explicit SimpleVacs(const FrmOnvField &mbf);

            const defs::inds_t &get();
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
            std::vector<defs::inds_t> m_inds;
            /**
             * map from the spinorb index to the "label"
             */
            const defs::inds_t m_map;
            /**
             * the simple decoding is cached at the same time, since it is often useful and available at small
             * additional cost
             */
            defs::inds_t m_simple_inds;

            LabelledBase(size_t nelement, const defs::inds_t &map, const FrmOnvField &mbf);

            const std::vector<defs::inds_t> &validated() const;

        public:
            /**
             * duplicate the spatial orbital irrep label map into two spin channels (spin major mapping)
             * @param grp_map
             *  spatial orbital irrep map
             * @return
             *  spin orbital irrep map with 2*nirrep labels
             */
            static defs::inds_t make_spinorb_map(const defs::inds_t &site_irreps, size_t nirrep);

            void clear();

            bool empty();

            size_t label(size_t ispinorb) const;
        };

        /**
         * update based on the set positions
         */
        struct LabelledOccs : LabelledBase {
        protected:
            LabelledOccs(size_t nelement, const defs::inds_t &map, const FrmOnvField &mbf);

        public:
            const std::vector<defs::inds_t> &get();

            const defs::inds_t &simple();
        };

        /**
         * update based on the clr positions
         */
        struct LabelledVacs : LabelledBase {
        protected:
            LabelledVacs(size_t nelement, const defs::inds_t &map, const FrmOnvField &mbf);

        public:
            const std::vector<defs::inds_t> &get();

            const defs::inds_t &simple();
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
            const std::vector<defs::inds_t> &m_inds_ref;

        public:
            NdBase(std::array<size_t, nind> shape, const std::vector<defs::inds_t> &inds) :
                    m_format(shape), m_inds_ref(inds) {}

            size_t size(const size_t &i) const {
                return m_inds_ref[i].size();
            }

            size_t size(const std::array<size_t, nind> &inds) const {
                return m_inds_ref[m_format.flatten(inds)].size();
            }

            const defs::inds_t &operator[](const size_t &i) const {
                return m_inds_ref[i];
            }

            const defs::inds_t &operator[](const std::array<size_t, nind> &inds) const {
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
            NdLabelledOccs(std::array<size_t, nind> shape, const defs::inds_t &map, const FrmOnvField &mbf) :
                    LabelledOccs(NdFormat<nind>(shape).m_nelement, map, mbf),
                    m_nd_inds(shape, m_inds) {}

            const NdBase<nind> &get() {
                LabelledOccs::get();
                return m_nd_inds;
            }

        protected:
            const NdBase<nind> &validated() const {
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
            NdLabelledVacs(std::array<size_t, nind> shape, const defs::inds_t &map, const FrmOnvField &mbf) :
                    LabelledVacs(NdFormat<nind>(shape).m_nelement, map, mbf),
                    m_nd_inds(shape, m_inds) {}

            const NdBase<nind> &get() {
                LabelledVacs::get();
                return m_nd_inds;
            }

        protected:
            const NdBase<nind> &validated() const {
                DEBUG_ASSERT_TRUE(is_valid(), "cache is not in sync with current MBF value");
                return m_nd_inds;
            }
        };

        /**
         * segregating the set or clear indices by spin value
         */
        struct SpinOccs : NdLabelledOccs<1> {
            SpinOccs(const FrmOnvField &mbf);
        };

        struct SpinVacs : NdLabelledVacs<1> {
            SpinVacs(const FrmOnvField &mbf);
        };

        /**
         * The immediate practical application of the general multidimensional decoding defined above, is in the case
         * where there is utility in segregating the set or clear indices by both spin and point group irrep.
         */
        struct SpinSymOccs : NdLabelledOccs<2> {
            SpinSymOccs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf);
        };

        struct SpinSymVacs : NdLabelledVacs<2> {
            SpinSymVacs(const AbelianGroupMap &grp_map, const FrmOnvField &mbf);
        };

        struct NonEmptyPairLabels : SimpleBase {
            NonEmptyPairLabels(const FrmOnvField &mbf);

            const defs::inds_t &get();
        };

        struct OccSites : SimpleBase {
            OccSites(const FrmOnvField &mbf);

            const defs::inds_t &get();
        };

        struct DoublyOccSites : SimpleBase {
            DoublyOccSites(const FrmOnvField &mbf);

            const defs::inds_t &get();
        };

        struct NotSinglyOccSites : SimpleBase {
            NotSinglyOccSites(const FrmOnvField &mbf);

            const defs::inds_t &get();
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
         * spin-partitioned occupied spin orbitals
         */
        frm::SpinOccs m_spin_occs;
        /**
         * spin-partitioned vacant spin orbitals
         */
        frm::SpinVacs m_spin_vacs;
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
        frm::OccSites m_occ_sites;
        /**
         * site indices with 2 electrons
         */
        frm::DoublyOccSites m_doubly_occ_sites;
        /**
         * site indices with 0 or 2 electrons
         */
        frm::NotSinglyOccSites m_not_singly_occ_sites;

        FrmOnv(const FrmOnvField &mbf);

        /**
         * clear all cached assets
         */
        void clear();
    };
}

#endif //M7_DECODEDFRMONV_H
