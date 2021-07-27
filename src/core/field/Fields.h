//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumberField.h"
#include "MultiField.h"
#include "FrmOnvField.h"
#include "BosOnvField.h"

namespace fields {
    template<typename T, size_t nind>
    using Numbers = NdNumberField<T, nind>;

    template<typename T>
    using Number = NumberField<T>;

    template<typename T, size_t nind>
    using NdBitset = BitsetField<T, nind>;
    template<typename T>
    using Bitset = NdBitset<T, 1ul>;

    template<size_t nind>
    using Flags = BitsetField<uint8_t, nind>;
    using Flag = BitField<uint8_t>;

    using FrmOnv = FrmOnvField;

    using BosOnv = BosOnvField;

    struct FrmBosOnv : MultiField<FrmOnv, BosOnv> {
        const std::string m_name;
        FrmOnv &m_frm;
        BosOnv &m_bos;

        FrmBosOnv(Row *row, size_t nsite, std::string name = "") :
                MultiField<FrmOnv, BosOnv>(row,
                                           {nullptr, nsite, name.empty() ? "" : name + " (fermion)"},
                                           {nullptr, nsite, name.empty() ? "" : name + " (boson)"}),
                m_name(name), m_frm(get<0>()), m_bos(get<1>()) {
        }

        FrmBosOnv(const FrmBosOnv &other) : FrmBosOnv(other.m_frm.row_of_copy(), other.m_frm.m_nsite,
                                                      other.m_name) {}

        FrmBosOnv &operator=(const FrmBosOnv &other) {
            m_frm = other.m_frm;
            m_bos = other.m_bos;
            return *this;
        }

        FrmBosOnv &operator=(const std::pair<defs::inds, defs::inds> &inds) {
            m_frm = inds.first;
            m_bos = inds.second;
            return *this;
        }
    };

//    struct FrmCsf : FrmOnv {
//        FrmCsf(Row *row, size_t nsite, std::string name = ""): FrmOnv(row, nsite, name){}
//    };

    typedef std::tuple<FrmOnv, FrmBosOnv, BosOnv> mbf_tup_t;

    template<size_t mbf_ind>
    using mbf_t = typename std::tuple_element<mbf_ind, mbf_tup_t>::type;
    typedef mbf_t<defs::mbf_ind> Mbf;

    struct FermionMevInds : MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> {
        typedef MultiField<Numbers<defs::mev_ind_t, 1>, Numbers<defs::mev_ind_t, 1>> base_t;
        const std::string m_name;
        Numbers<defs::mev_ind_t, 1> &m_ann;
        Numbers<defs::mev_ind_t, 1> &m_cre;

        FermionMevInds(Row *row, size_t nann, size_t ncre, std::string name = "") :
                base_t(row, {nullptr, {nann}, name.empty() ? "" : name + "_ann"},
                       {nullptr, {ncre}, name.empty() ? "" : name + "_cre"}),
                m_name(name), m_ann(get<0>()), m_cre(get<1>()) {
        }

        FermionMevInds(const FermionMevInds &other) :
                FermionMevInds(other.m_ann.row_of_copy(),
                               other.m_ann.m_format.m_nelement,
                               other.m_cre.m_format.m_nelement, other.m_name) {}

        FermionMevInds &operator=(const FermionMevInds &other) {
            m_ann = other.m_ann;
            m_cre = other.m_cre;
            return *this;
        }

        FermionMevInds &operator=(const std::pair<defs::inds, defs::inds> &inds) {
            for (size_t i = 0ul; i < inds.first.size(); ++i) m_ann[i] = inds.first[i];
            for (size_t i = 0ul; i < inds.second.size(); ++i) m_cre[i] = inds.second[i];
            return *this;
        }

        /**
         * The RDM entry keys represent the ascending-ordered sector of the normal-ordered pure expectation values,
         * the full spin-resolved RDM can be constructed from this via permutations of the indices.
         * @return
         *  true if the SQ creation and annihilation operator index vectors are properly ordered
         */
        bool is_ordered() const {
            return m_ann.is_ordered(false, true) && m_cre.is_ordered(false, true);
        }

        /**
         * all elements of the RDM have the same rank, but not the same excitation level. If the number of indices in
         * common between the creation and annihilation operators is zero, the excitation level is the same as the rank
         * @return
         *  number of indices in common between ascending ordered creation and annihilation spin orbital operator strings
         */
        size_t ncommon_sq_op_ind() const {
            size_t ncommon = 0ul;
            size_t icre = 0ul;
            size_t iann = 0ul;
            ASSERT(is_ordered());
            while(icre<m_cre.nelement() && iann<m_ann.nelement()){
                if (m_cre[icre]>m_ann[iann]) ++icre;
                else if (m_cre[icre]<m_ann[iann]) ++iann;
                else {
                    // common element found
                    ++ncommon;
                    ++icre; ++iann;
                }
            }
            return ncommon;
        }

        void common_sq_op_inds(defs::inds& common) const {
            common.clear();
            size_t icre = 0ul;
            size_t iann = 0ul;
            ASSERT(is_ordered());
            while(icre<m_cre.nelement() && iann<m_ann.nelement()){
                if (m_cre[icre]>m_ann[iann]) ++icre;
                else if (m_cre[icre]<m_ann[iann]) ++iann;
                else {
                    // common element found
                    common.push_back(m_cre[icre]);
                    ++icre; ++iann;
                }
            }
        }
    };

    /**
     * Fields are defined as symbols, i.e. members within a Row-derived class.
     * They also retain a pointer to the Row with which they are associated, and the Row in turn keeps a vector of
     * corresponding pointers to FieldBase objects. In some situations, we may want to specify a field within a row,
     * and them use that field but in the context of another row. In order to do this, the offset in memory between
     * the specified row and field is computed, then that offset is applied to the pointer to the other row, resulting
     * in a reference to the corresponding field in the other row.
     * This function is subject to run time checks ensuring safety of the casting operations.
     * @tparam row_t
     *  Row-derived class containing the field_t symbol
     * @tparam field_t
     *  FieldBase-derived class defined as a member of row_t
     * @param target
     *  Row object for which the field reference is desired
     * @param source
     *  Row object for which the field reference is specified
     * @param field
     *  FieldBase-derived object which is the referenced symbol within the source row
     * @return
     *  reference to the same field within "target" as is represented by "field" within "source"
     */
    template<typename row_t, typename field_t>
    static field_t& identify(row_t& target, row_t& source, field_t& field){
        static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
        static_assert(std::is_base_of<FieldBase, field_t>::value, "Template arg must be derived from FieldBase");
        REQUIRE_TRUE(static_cast<FieldBase&>(field).belongs_to_row(source), "field arg must belong to source arg");
        auto target_ptr = reinterpret_cast<char*>(&target);
        auto source_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&field);
        long byte_offset = field_ptr - source_ptr;
        DEBUG_ASSERT_GT(byte_offset, 0l, "field pointer is not positively offset from row!");
        auto ptr = reinterpret_cast<field_t*>(target_ptr + byte_offset);
        DEBUG_ASSERT_TRUE(ptr->belongs_to_row(target), "field identification failed");
        return *ptr;
    }

    /**
     * same as the above function, but for references to MultiFields instead
     * @tparam row_t
     *  Row-derived class containing the field_t symbol
     * @tparam Args
     *  FieldBase-derived arguments to the MultiField class template
     * @param target
     *  Row object for which the multifield reference is desired
     * @param source
     *  Row object for which the multifield reference is specified
     * @param multifield
     *  MultiField-derived object which is the referenced symbol within the source row
     * @return
     *  reference to the same MultiField within "target" as is represented by "field" within "source"
     */
    template<typename row_t, typename ...Args>
    static MultiField<Args...>& identify(row_t& target, row_t& source, MultiField<Args...>& multifield){
        static_assert(std::is_base_of<Row, row_t>::value, "Template arg must be derived from Row");
        REQUIRE_TRUE(multifield.belongs_to_row(source), "multifield arg must belong to source arg");
        auto target_ptr = reinterpret_cast<char*>(&target);
        auto source_ptr = reinterpret_cast<char*>(&source);
        auto field_ptr = reinterpret_cast<char*>(&multifield);
        long byte_offset = field_ptr - source_ptr;
        DEBUG_ASSERT_GT(byte_offset, 0, "field pointer is not positively offset from row!");
        auto ptr = reinterpret_cast<MultiField<Args...>*>(target_ptr + byte_offset);
        DEBUG_ASSERT_TRUE(ptr->belongs_to_row(target), "multifield identification failed");
        return *ptr;
    }

}

#endif //M7_FIELDS_H